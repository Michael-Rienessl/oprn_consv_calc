'''
@author: ichaudr

'''

from Bio import Entrez, SeqIO
from features import GenomeFeature
from operon import Operon
import time
import os
import logging

class GenomeFragment:
    '''
    Represents the hits associated with a specific genome fragment (i.e. plasmind, chromosome, contig). The representative operon(s) can be assembled from the set of hits associated with the fragment.
    
    '''
    
    def __init__(
        self,
        name,
        genome_fragment_accession,
        req_limit,
        sleep_time,
        cache_directory,
        allow_ncbi_download=True
    ):
        
        self.cache_directory = cache_directory
        self.allow_ncbi_download = allow_ncbi_download

        # Optional direct local GenBank file path.
        # Set by operon_conserve_detect.py in local BLAST mode.
        self.local_gbk_file_path = None

        self.hits = []
        self.all_features = []
        self.operons = []

        self.name = name
        self.genome_accession = genome_fragment_accession
        self.assembly_accession = None
        self.taxid = 'None'
        self.req_limit = req_limit
        self.sleep_time = sleep_time

        self.full_record = None
        self.species_name = None

    def fetch_features(self):
        '''
        Obtain all coding features in the genome and save them as GenomeFeature objects.
        Uses a SeqRecord-based parser for both local and remote records.
        '''

        if self.full_record is None:
            self.fetch_record()
        
        if self.full_record is None:
            self.all_features = []
            return

        for feat in self.full_record.features:
            gf = self._seqfeature_to_genomefeature(feat, self.full_record)
            if gf is not None:
                self.all_features.append(gf)

        self.sort_all_features()

        
        logging.info("=== GENBANK PARSER SANITY CHECK ===")
        logging.info(f"Fragment Accession: {getattr(self, 'genome_accession', 'N/A')}")
        logging.info(f"Total features parsed: {len(self.all_features)}")


    def sort_all_features(self):
        '''
        Sorts the all_features list in-place from 5' to 3'.
        Should be called once after all features have been loaded.
        '''
        self.all_features.sort(key=lambda x: x.five_end)

    def fetch_record(self):
        """
        Obtains the full GenBank record for the fragment.

        Local mode:
            - uses self.local_gbk_file_path if provided,
            - never downloads from NCBI if allow_ncbi_download=False,
            - does not modify the accession,
            - expects self.genome_accession == GenBank record.id.

        Remote mode:
            - uses the fast cache convention first:
                  <cache_directory>/<genome_accession>.gb
            - if this exact cache file exists, it is opened directly.
            - if it does not exist, the record is downloaded from NCBI and
              written to that exact file name.
        """
        import os
        import logging
        from Bio import Entrez, SeqIO

        target_acc = str(self.genome_accession).strip()

        # ------------------------------------------------------------
        # 1. Fast local path from pre-built index
        # ------------------------------------------------------------
        local_path = getattr(self, "local_gbk_file_path", None)

        if local_path and os.path.exists(local_path):
            try:
                for record in SeqIO.parse(local_path, "genbank"):
                    record_id = str(record.id).strip()

                    if record_id == target_acc:
                        self.full_record = record
                        self.species_name = record.annotations.get("organism", "Unknown Species")

                        logging.info(
                            f"Loaded indexed local record for {self.genome_accession} "
                            f"from {os.path.basename(local_path)}"
                        )
                        return

                logging.error(
                    f"Indexed local file exists, but no exact record.id match was found. "
                    f"target={target_acc}, file={local_path}"
                )

                self.full_record = None
                self.species_name = "Unknown Species"
                return

            except Exception as e:
                logging.error(f"Failed to read indexed local GenBank file {local_path}: {e}")
                self.full_record = None
                self.species_name = "Unknown Species"
                return

        # ------------------------------------------------------------
        # 2. Local mode safety lock
        # ------------------------------------------------------------
        if not self.allow_ncbi_download:
            logging.error(
                f"[LOCAL MODE BLOCKED] No indexed local GBK file for {self.genome_accession}. "
                f"NCBI download is disabled."
            )
            self.full_record = None
            self.species_name = "Unknown Species"
            return

        # ------------------------------------------------------------
        # 3. Remote mode: direct cache lookup by accession-based filename
        # ------------------------------------------------------------
        os.makedirs(self.cache_directory, exist_ok=True)

        # Remote cache convention:
        #     Genome accession CP000000.1 -> cache/CP000000.1.gb
        record_file = os.path.join(self.cache_directory, f"{target_acc}.gb")

        if os.path.exists(record_file):
            try:
                with open(record_file, "r", encoding="utf-8") as handle:
                    self.full_record = SeqIO.read(handle, "genbank")

                self.species_name = self.full_record.annotations.get("organism", "Unknown Species")
                logging.info(
                    f"Loaded remote cached record for {self.genome_accession} "
                    f"from {os.path.basename(record_file)}"
                )
                return

            except Exception as e:
                logging.warning(
                    f"Cached GenBank file exists but could not be parsed: "
                    f"{record_file}. Will attempt to re-download. Error: {e}"
                )
                self.full_record = None

        # Compatibility direct checks for cache files written by older/newer
        # pipeline versions. These are still O(1)-style checks because they
        # test exact expected filenames, not the whole cache directory.
        for ext in (".gbk", ".gbff", ".genbank"):
            alt_file = os.path.join(self.cache_directory, f"{target_acc}{ext}")

            if not os.path.exists(alt_file):
                continue

            try:
                with open(alt_file, "r", encoding="utf-8") as handle:
                    self.full_record = SeqIO.read(handle, "genbank")

                self.species_name = self.full_record.annotations.get("organism", "Unknown Species")
                logging.info(
                    f"Loaded remote cached record for {self.genome_accession} "
                    f"from {os.path.basename(alt_file)}"
                )
                return

            except Exception as e:
                logging.warning(
                    f"Alternative cached GenBank file exists but could not be parsed: "
                    f"{alt_file}. Error: {e}"
                )
                self.full_record = None

        # ------------------------------------------------------------
        # 4. Remote NCBI fallback: download and save as <accession>.gb
        # ------------------------------------------------------------
        logging.info(f"Downloading GenBank record for {self.genome_accession} from NCBI...")

        record_text = None

        for attempt in range(self.req_limit):
            try:
                handle = Entrez.efetch(
                    db="nuccore",
                    id=target_acc,
                    strand=1,
                    seq_start="begin",
                    seq_stop="end",
                    rettype="gbwithparts",
                    retmode="text"
                )

                record_text = handle.read()
                handle.close()
                time.sleep(self.sleep_time)
                break

            except Exception as e:
                logging.warning(
                    f"NCBI exception on attempt {attempt + 1}/{self.req_limit} "
                    f"for {self.genome_accession}: {e}"
                )
                time.sleep(self.sleep_time)

        if not record_text:
            logging.error(
                f"Could not download GenBank record for {self.genome_accession} "
                f"after {self.req_limit} attempts."
            )
            self.full_record = None
            self.species_name = "Unknown Species"
            return

        try:
            with open(record_file, "w", encoding="utf-8") as handle:
                handle.write(record_text)

            with open(record_file, "r", encoding="utf-8") as handle:
                self.full_record = SeqIO.read(handle, "genbank")

            self.species_name = self.full_record.annotations.get("organism", "Unknown Species")
            logging.info(
                f"Downloaded and cached remote record for {self.genome_accession} "
                f"as {os.path.basename(record_file)}"
            )

        except Exception as e:
            logging.error(
                f"Error writing or parsing downloaded GenBank record for "
                f"{self.genome_accession}: {e}"
            )
            self.full_record = None
            self.species_name = "Unknown Species"

    def fetch_hit_features(self, margin_limit=20, max_attempts=5, mult_factor=3):
        '''
        Attempts to fetch the full feature for each of the hits assigned to this genome fragment. Filters the list of hits so that there is no more than 1 hit per feature.

        Parameters
        ----------
        None

        Returns
        -------
        None - features are set internally for each hit object. 
        '''

        # If file doesn't exist in local mode we don't pass None, but abbort.
        if self.full_record is None:
            import logging
            logging.error(f"Skipping feature mapping for {self.genome_accession} - Record not available locally.")
            return

        #Fetch the feature for each hit
        for hit in self.hits:
            hit.fetch_feature(self.full_record)
            #hit.fetch_feature(self.full_record, margin_limit=margin_limit, max_attempts=max_attempts, mult_factor=mult_factor)    

    def _seqfeature_to_genomefeature(self, feat, record):
        """
        Convert a Biopython SeqFeature (CDS) into a GenomeFeature object.
        Returns None for non-CDS or unprocessable features.
        """

        if feat.type != 'CDS':
            return None
        
        qualifiers = feat.qualifiers if feat.qualifiers else {}

        locus_tag = qualifiers.get('locus_tag', ['None'])[0]
        protein_accession = qualifiers.get('protein_id', ['None'])[0]

        if protein_accession in [None, 'None', ''] and locus_tag not in [None, 'None', '']:
            protein_accession = locus_tag

        # Prefer annotated translation, otherwise try translating from nucleotide.
        aa_sequence = qualifiers.get('translation', [None])[0]
        if not aa_sequence:
            try:
                aa_sequence = str(feat.extract(record.seq).translate(to_stop=True))
            except Exception:
                aa_sequence = 'None'

        coding_start = int(feat.location.start) + 1
        coding_end = int(feat.location.end)

        strand_val = feat.location.strand
        if strand_val == 1:
            strand = '+'
            five_end = coding_start
            three_end = coding_end
        elif strand_val == -1:
            strand = '-'
            five_end = coding_end
            three_end = coding_start
        else:
            strand = '?'
            five_end = coding_start
            three_end = coding_end
        
        if protein_accession == None:
            protein_accession = "None"
        
        if locus_tag == None:
            locus_tag = "None"

        if aa_sequence == None:
            aa_sequence = "None"

        gf = GenomeFeature(
            genome_accession=self.genome_accession,
            genome_fragment_name=self.name,
            req_limit=self.req_limit,
            sleep_time=self.sleep_time,
            strand=strand,
            aa_sequence=aa_sequence,
            coding_start=coding_start,
            coding_end=coding_end,
            five_end=five_end,
            three_end=three_end,
            protein_accession=protein_accession,
            locus_tag=locus_tag
        )

        # Optional annotation
        if 'product' in qualifiers and qualifiers['product']:
            gf.product = qualifiers['product'][0]
        
        if 'gene' in qualifiers and qualifiers['gene']:
            gf.gene = qualifiers['gene'][0]
        
        return gf


    def purge_hits(self):
        '''
        Removes any duplicates from the list of hits.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        
        #Remove any duplicate features so that no more than
        logging.info(f"Purging hits... Pre-purge count: {len(self.hits)}")  

        seen_accessions = set()
        purged_hits = []
        for hit in self.hits:
            if hit.protein_accession not in seen_accessions:
                purged_hits.append(hit)
                seen_accessions.add(hit.protein_accession)
            else:
                logging.debug(f"DUPLICATE REMOVED: {hit.protein_accession} on {hit.genome_accession} (Locus Tag: {hit.locus_tag})")
        
        self.hits = purged_hits
        logging.info(f"Purging hits... Post-purge count: {len(self.hits)}")
            
        
    
    ####NOT USED#####
    def fetch_taxid(self):
        '''
        Fetches the taxonomic id for this fragment

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        
        #Fetches a much smaller region of the full record to parse out the taxid
        for i in range(self.req_limit):

            try:

                handle = Entrez.efetch(db="nuccore", id=self.genome_accession, seq_start=0, seq_stop=1, rettype='gb', retmode='XML')
                record = Entrez.read(handle, 'xml')
                time.sleep( self.sleep_time)

                break

            except:

                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                if i == (self.req_limit - 1):
                        print("\t\tCould not download record after " + str(self.req_limit) + " attempts")
        
        for feature in record[0]['GBSeq_feature-table']:
            #Get taxid
            if feature['GBFeature_key'] == 'source':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'db_xref':
                        self.taxid = qual['GBQualifier_value'].split(':')[1]
                    if qual['GBQualifier_name'] == 'organism':
                        self.species_name = qual['GBQualifier_value']
        
        #Clear up memory
        del record
        record = None

    def add_hit(self, a_hit):
        '''
        Adds an AnnotatedHit object to the list of hits that belong to this genome fragment.

        Parameters
        ----------
        a_hit : AnnotatedHit object
            A hit that is associated with this genome fragment and will be added to the fragment's record.

        Returns
        -------
        None.

        '''
        
        #Check if the AnnotatedHit that is passed belongs with this fragment
        if a_hit.genome_fragment_name == self.name:
            self.hits.append(a_hit)

    def assemble_operons(self, feature_limit=3, intergenic_limit=150):
        '''
        Takes the list of hits and organizes them into putative operons by following the scheme:
            1. Separate hits that are on + and - strands. 
            2. Sort the hits from 5' to 3' 
            3. Fill in gaps with features only if they are on the same strand within some intergenic distance, intergenic_limit. A maximum of feature_limit features will be added to fill the gaps.                            
            4. After filling in all the gaps, set of extended hits are either grouped as a single operon or split depending on the intergenic distance limit, intergenic_limit. 
        
        Parameters
        ----------
        feature_limit: int
            The max number of features that are allowed to space two genes.
        intergenic_limit: int
            The max distance allowed between two genes.
        
        Returns
        -------
        None - operons are added the fragment list
        '''
        
        #Hold the hits on the minus and plus strands
        plus_strand_hits = []
        minus_strand_hits = []

        #Group all hits by strand
        for hit in self.hits:
            if hit.strand == '+':
                plus_strand_hits.append(hit)
            elif hit.strand == '-':
                minus_strand_hits.append(hit)
        
        #Sort plus and minus strand hits from 5' to 3'
        if len(plus_strand_hits) > 0:
            plus_strand_hits = sorted(plus_strand_hits, key=lambda hit: hit.five_end)
        
        if len(minus_strand_hits) > 0:
            minus_strand_hits = sorted(minus_strand_hits, key=lambda hit: hit.five_end)


        #Start assembling the operons on each of the strands

        while len(plus_strand_hits) > 0:
            temp_operon_features = []
            temp_operon_features.append(plus_strand_hits.pop(0))

            #Keep track of inserted intergenic features added
            num_features_added = 0

            continue_adding = True
            while continue_adding:             

                #Get the next feature in the genome
                next_feat = self.get_next_feature(temp_operon_features[-1])

                #Keep track of whether or not the next hit has been added
                added_next_feat = False

                if next_feat == None:
                    continue_adding = False
                    continue

                #Check if strand is valid
                if next_feat.strand == '+':

                    #If the next feature is equal to the next hit in the plus strand list, append it and reset the reference and number of features added
                    if len(plus_strand_hits) > 0:
                        if next_feat == plus_strand_hits[0]:
                            temp_operon_features.append(plus_strand_hits.pop(0))
                            num_features_added = 0
                            added_next_feat = True
                    
                    #If the next feature is not in the original list of hits (i.e. its an inserted intergenic element), check if it meets the intergenic disantce
                    # and add to the max number of intergenic features. 
                    if not added_next_feat:
                        if temp_operon_features[-1].get_intergenic_distance(next_feat) <= intergenic_limit and num_features_added <= feature_limit:
                            temp_operon_features.append(next_feat)
                            num_features_added = num_features_added + 1
                            added_next_feat = True
                        else:
                            continue_adding = False


                #If strand is not valid, break the operon    
                else:
                    continue_adding = False
            
            self.add_operon(temp_operon_features, "+")

        while len(minus_strand_hits) > 0:
            temp_operon_features = []
            temp_operon_features.append(minus_strand_hits.pop(0))

            #Keep track of inserted intergenic features added
            num_features_added = 0

            continue_adding = True
            while continue_adding == True:

                #Get the next feature in the genome
                next_feat = self.get_next_feature(temp_operon_features[-1])

                #Keep track of whether or not the next hit has been added
                added_next_feat = False

                if next_feat == None:
                    continue_adding = False
                    continue

                #Check if strand is valid
                if next_feat.strand == '-':

                    #If the next feature is equal to the next hit in the plus strand list, append it and reset the reference and number of features added
                    if len(minus_strand_hits) > 0:
                        if next_feat == minus_strand_hits[0]:
                            temp_operon_features.append(minus_strand_hits.pop(0))
                            num_features_added = 0
                            added_next_feat = True
                    
                    #If the next feature is not in the original list of hits (i.e. its an inserted intergenic element), check if it meets the intergenic disantce
                    # and add to the max number of intergenic features. 
                    if not added_next_feat:
                        if temp_operon_features[-1].get_intergenic_distance(next_feat) <= intergenic_limit and num_features_added <= feature_limit:
                            temp_operon_features.append(next_feat)
                            num_features_added = num_features_added + 1
                            added_next_feat = True
                        else:
                            continue_adding = False


                #If strand is not valid, break the operon    
                else:
                    continue_adding = False
            
            self.add_operon(temp_operon_features, "-")

    def get_next_feature(self, current_feature):
        '''
        Get the next feature in the genome that is on the same strand (always considers it from 5'-3' regardless of strand)

        Parameters
        ----------
        current_feature: GenomeFeature object
            The feature to move forward from. 
        
        Returns
        -------
        next_feature: GenomeFeature object
            The next sequential feature 5' to 3'
        '''

        #Get the position of the current feature
        pos = self.get_feature_position(current_feature)

        if pos == -1:
            print("Feature: " + str(current_feature) + "\nNot found in " + str(self.name))
            return None

        if pos < len(self.all_features) - 1:
            return self.all_features[pos + 1]
        else:
            return None

    def get_feature_position(self, feature):
        '''
        Gets the position of a feature in the genome_features list

        Parameters
        ----------
        feature: GenomeFeature object
            Feature of interest
        
        Returns
        -------
        positon: int
            Index of feature of interest, -1 if not found
        '''

        position = 0
        found = False

        while position < (len(self.all_features)):
            if feature == self.all_features[position]:
                found = True
                break

            position = position + 1
        
        if found == True:
            return position
        else:
            return -1

    def add_operon(self, feature_list, strand):
        '''
        Appends a new operon object to the genome fragment.

        Parameters
        ----------
        feature_list: list[GenomeFeature object]
            List of features belonging to the operon
        strand: string
            The strand, + or -, that the operon is located on.
        
        Returns
        -------
        None
        '''
        self.operons.append(Operon(genome_fragment_name = self.name, genome_accession=self.genome_accession ,genome_features=self.all_features, strand=strand))

        for f in feature_list:
            self.operons[-1].add_feature(f)
    
    ####NOT USED#####
    def get_species_name(self, parse_phrases=['whole', 'plasmid', 'chromosome', 'contig', 'NODE', 'scaffold', 'complete', 'Contig', 'Scaffold', 'sequence', 'Sequence', ]):
        '''
        Parses the species name from the GenomeFrament name that comes from the BLAST hit discription

        Paramters
        ---------
        parse_phrases: list[string]
            The different phrases to cutoff at in the hit description to parse out the species name

        Returns
        -------
        None
        '''

        if len(self.species_name) > 0:
            return
        
        ###Testing
        for i in range(self.req_limit):

            try:

                handle = Entrez.efetch(db="taxonomy", id=self.taxid)
                record = handle.read()
                time.sleep(self.sleep_time)

                break

            except:
                print('Messed up here')
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                if i == (self.req_limit - 1):
                        print("\t\tCould not download record after " + str(self.req_limit) + " attempts")
        
        self.species_name = record[(record.find('<ScientificName>') + 16):record.find('</ScientificName>')]

        del record
        record = None

        '''
        temp_name = self.name

        #Holds the locations where each phrase appears in the fragment name
        phrase_locations = []

        #Get locations of each phrase
        for phrase in parse_phrases:

            #Location of the phrase
            phrase_pos = temp_name.find(phrase)

            if phrase_pos != -1:
                phrase_locations.append((phrase, phrase_pos))
        
        #Split the fragment name at the first phrase that occurs in the name 
        if len(phrase_locations) > 0:

            #Stores the split position and the phrase to split at
            split_pos = phrase_locations[0][1]
            split_phrase = ""
            
            for pair in phrase_locations:
                if pair[1] <= split_pos:
                    split_pos = pair[1]
                    split_phrase = pair[0]
            if split_pos != -1:
                self.species_name = temp_name.split(split_phrase)[0][0:-1] '''           

    def __str__(self):
        to_return = "Genome Fragment Name: " + self.name + "\n" + "Total number of hits: " + str(len(self.hits)) + "\nTotal number of features: " + str(len(self.all_features))  + "\nTotal number of operons: " + str(len(self.operons)) + "\n"
        
        for operon in self.operons:
            to_return = to_return + str(operon) + "\n"

        return to_return    