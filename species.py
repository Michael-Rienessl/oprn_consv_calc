'''
@author: ichaudr

'''

from Bio.Graphics.GenomeDiagram import Diagram, FeatureSet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from features import AnnotatedHit
import re
import os
import math
import logging

class Species:
    '''
    Holds all the GenomeFragment objects that belong to the same species. 
    '''

    def __init__(self, assembly_accession, genome_fragments=[]):
        self.assembly_accession = assembly_accession
        self.species_name = 'No Name Assigned'
        self.genome_fragments = genome_fragments
        self.genome_fragments_accessions = []
        self.sim_score = 0
        self.query_percent_ids = {}
        self.taxid = '-1'
    
    def calculate_total_avg_aai(self):
        """
        Calculates the average AAI of ALL hits found in the genome for this species,
        not just the average of the maximums.
        """
        all_identities = []
        for frag in self.genome_fragments:
            for hit in frag.hits:
                if hasattr(hit, 'percent_identity'):
                    all_identities.append(float(hit.percent_identity))
        
        if all_identities:
            return sum(all_identities) / len(all_identities)
        return 0.0
    


    ###This function is not used anywhere anymore
    def contains_query_feature(self, query_accession):
        '''
        Searches all the GenomeFragment objects for a feature with that resulted from an original query from the reference operon

        Parameters
        ----------
        query_accession: string
            The accession of the query of interest
        
        Returns
        -------
        found_feature: GenomeFeature object
            The GenomeFeature object that resutled from the query of interest
        '''

        found_feature = None

        for fragment in self.genome_fragments:
            for hit in fragment.hits:
                if hit.query_accession == query_accession:
                    found_feature = hit
        
        return found_feature
    

    def get_query_percent_ids(self, query_accessions):
        '''
        Determines the hit with the max percent identity (per BLAST search) for a given list query accessions and stores it in a dictionary.

        Parameters
        ----------
        query: list[string]
            The accessions for the query from the reference operon.
        
        Returns
        -------
        None - a dictionary is set in the object with {key:value} = {reference_accession:percent_id}
        '''
        
        for q_acc in query_accessions:
            #Stores the max percenty identity 
            max_percent_id = 0

            #Iterate through all genome fragments and all features to find any corresponding with the target query and assess its percent identity. 
            for frag in self.genome_fragments:
                for hit in frag.hits:
                    if hit.query_accession == q_acc:
                        if hit.percent_identity > max_percent_id:
                            max_percent_id = hit.percent_identity
            
            self.query_percent_ids[q_acc] = max_percent_id


    
    def add_genome_fragment(self, genome_fragment):
        '''
        Adds a GenomeFragment to the Species object.

        Parameters
        ----------
        genome_fragment: GenomeFragment object
            The GenomeFragment to append to the genome_fragments list for this Species object
        
        Returns
        -------
        None
        '''

        if genome_fragment.assembly_accession == self.assembly_accession:
            self.genome_fragments.append(genome_fragment)


    def measure_sim(self, reference_operon):
        '''
        Calculates two types of structural conservation:
        1. Total Genomic SIM: A global score normalized by gene frequency (Paper Logic).
        2. Best Operon SIM: Identifies the single best-preserved operon fragment.

        Detailed Workflow:
            1. Pair Generation: Creates unique pairs from the reference_operon accessions.
            2. Global Gene Counting: Counts every occurrence of each query gene across 
               the entire genome to establish the "Max Possible" denominator.
            3. Fragment Evaluation: Iterates through all detected operon fragments:
                - Calculates "Local SIM" and "Local AAI" for each fragment.
                - Identifies the "Best Operon" (highest Local SIM, then highest Local AAI).
            4. Pair Tallies: Collects all adjacent gene pairs found in all fragments.
            5. Total SIM Calculation: For each reference pair, calculates:
               (Found instances in all fragments) / min(Count of Gene A, Count of Gene B).
            6. Final Scoring: Averages these weighted values by the number of reference pairs.

        Parameters
        ----------
        reference_operon: list[string]
            Accessions of the query proteins in the exact reference order.

        Returns
        -------
        float: The Total Genomic SIM score (also stored in self.sim_score).
        '''

        # 1. Setup Reference Pairs
        ref_pairs = []
        for i in range(len(reference_operon)-1):
            ref_pairs.append((reference_operon[i], reference_operon[i+1]))

        logging.info(f"    [Struct] Perfect Reference Pairs: {ref_pairs}")

        # 2. Global Gene Counting (for Total Genomic SIM normalization)
        species_query_count = {acc: 0 for acc in reference_operon}
        for frag in self.genome_fragments:
            for hit in frag.hits:
                if hit.query_accession in species_query_count:
                    species_query_count[hit.query_accession] += 1
        
        # Store for summarized output
        self.query_hits_counts = species_query_count

        # 3. Analyze Fragments and Find Best Operon
        species_pairs = []
        self.best_operon = None
        max_local_sim = -1.0
        max_local_aai = -1.0

        for frag in self.genome_fragments:
            for operon in frag.operons:
                # Use the new method from operon.py (Step 1)
                operon.calculate_local_stats(ref_pairs)
                
                # Logic to identify the Best Operon in the species
                if (operon.local_sim > max_local_sim) or \
                   (operon.local_sim == max_local_sim and operon.local_aai > max_local_aai):
                    max_local_sim = operon.local_sim
                    max_local_aai = operon.local_aai
                    self.best_operon = operon

                # Collect observed pairs for the Total SIM math
                # We only count pairs consisting of reference genes
                purged = [f for f in operon.features if hasattr(f, 'query_accession')]
                if len(purged) > 1:
                    for i in range(len(purged) - 1):
                        species_pairs.append((purged[i].query_accession, purged[i+1].query_accession))

        # 4. Total Genomic SIM Calculation (Paper-based)
        weighted_num_match = 0
        for ref_pair in ref_pairs:
            # Count how many times this pair appears (both directions if needed)
            num_occur = species_pairs.count(ref_pair) + species_pairs.count((ref_pair[1], ref_pair[0]))
            
            # Max possible is the "bottleneck" (minimum count of the two genes)
            max_possible_occur = min(species_query_count[ref_pair[0]], species_query_count[ref_pair[1]])

            if max_possible_occur > 0:
                added_weight = (num_occur / max_possible_occur)
                weighted_num_match += added_weight
                logging.info(f"    [Struct] Pair {ref_pair} Math: Found {num_occur} / Max {max_possible_occur} -> Weight: {added_weight:.2f}")
            else:
                logging.info(f"    [Struct] Pair {ref_pair} Math: Max possible is 0 -> Weight: 0.00")

        # 5. Finalize
        self.sim_score = weighted_num_match / len(ref_pairs) if ref_pairs else 0.0
        return self.sim_score

    def extract_taxids(self):
        '''
        Gets the tax ID for the species by extracting the largest taxid from the fragments.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''

        current_taxid = int(self.genome_fragments[0].taxid)

        for frag in self.genome_fragments:
            frag_taxid = int(frag.taxid)

            if frag_taxid > current_taxid:
                current_taxid = frag_taxid
        
        self.taxid = current_taxid
    
    def extract_genome_accessions(self):
        '''
        Gets the genome accessions for all the genome fragments.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        accessions_raw = []
        for frag in self.genome_fragments:
            accessions_raw.append(frag.genome_accession)
        
        accessions_purged = []
        for a in accessions_raw:
            if not(a in accessions_purged):
                accessions_purged.append(a)
        
        self.genome_fragments_accessions = accessions_purged

    def extract_names(self):
        '''
        Extracts the species name - picks the shortest name from the list of fragments.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        
        current_name = self.genome_fragments[0].species_name

        for frag in self.genome_fragments:
            frag_name = frag.species_name

            if frag_name is None:
                continue

            if current_name is None or current_name == "Unknown Species":
                current_name = frag_name

            elif len(frag_name) < len(current_name):
                current_name = frag_name

        if current_name is None:
            self.species_name = f"Unknown_Species_{self.assembly_accession}"
        else:
            self.species_name =  re.sub("[^0-9a-zA-Z]+", "_", current_name)
    

    def draw_figure(self, color_code={}, output_dir='./output/svg/'):
        '''
        Draws a SVG map of the operons that were detected. 

        Parameters
        ----------
        color_code: {reference_gene:color}
            The colors to be used for hits corresponding to each of the reference genes. 

        Returns
        -------
        None
        '''
        
        #Make the output directory for this species
        output_dir = output_dir + self.species_name + '/'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        for fragment in self.genome_fragments:

            #Keep track of the operon count to name the tracks
            operon_count = 0
            
            #Iterate through the features in each of the operons in the current fragment and add them to the current feature set.
            for operon in fragment.operons:

                #The GenomeDiagram for this operon
                diagram = Diagram()

                #Increment the operon_count
                operon_count = operon_count + 1

                #Keep track of the start/end of the track
                if len(operon.features) == 0:
                    continue
                t_start = int(operon.features[0].five_end)
                t_end = int(operon.features[0].three_end)

                current_track = diagram.new_track(1, greytrack=1, greytrack_labels=2, greytrack_fontsize=3, name=(fragment.genome_accession+ '|\noperon' +str(operon_count)), scale=1, scale_ticks=0)

                #Make a feature set for this track
                curr_feat_set = current_track.new_set()


                for feat in operon.features:

                    #Pull all the paramaters for the current feature
                    raw_start = int(feat.five_end)
                    raw_end = int(feat.three_end)
                    
                    # BIOPYTHON FIX: start must ALWAYS be <= end for FeatureLocation
                    start = min(raw_start, raw_end)
                    end = max(raw_start, raw_end)

                    # Adjust the start/end of the diagram if needed
                    if start < t_start:
                        t_start = start
                    if end > t_end:
                        t_end = end


                    strand = 0
                    if feat.strand == '+':
                        strand = 1
                    elif feat.strand == '-':
                        strand = -1


                    gene_name = feat.protein_accession

                    #Determine the color of the feature
                    color = ''
                    if not isinstance(feat, AnnotatedHit):
                        color = color_code['intergenic']
                    else:
                        color = color_code[feat.query_accession]

                    #Make the SeqFeature object
                    curr_feat = SeqFeature(location=FeatureLocation(start=start, end=end, strand=strand), qualifiers={'gene':[gene_name]})

                    #Add the current feature to the current feature set
                    curr_feat_set.add_feature(curr_feat, label=True, label_size=4, color=color, sigil='ARROW')


                diagram.draw(format='linear', pagesize=(150,100), fragments=1, tracklines=False, fragment_size=.25, start=int(t_start)-20, end=(t_end)+20, xr=0, xl=0, yt=0, yb=0)

                #Make the filename
                filename = str(output_dir + fragment.genome_accession + '_' + str(operon_count) + ".svg")

                #Write the output
                diagram.write(filename=filename, output='SVG')


    def clean(self):
        '''
        Deletes the genome fragments stored in this object to free up memory. 
        '''
        del self.genome_fragments
        self.genome_fragments = None


