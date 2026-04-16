'''
@author: ichaudr

Determines how conserved a specific reference operon is within a specific taxonomic clade. 

'''

from Bio import Entrez, SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from features import AnnotatedHit, GenomeFeature
from genome_fragment import GenomeFragment
from species import Species
from tqdm import tqdm
import datetime
import time
import json
import csv
import glob
import re
import sys
import csv
import threading
import uuid
import os
import shutil
import subprocess
import logging
import os
import socket

socket.setdefaulttimeout(60)


# --- LOGGING SETUP ---
# Erstellt eine Datei "pipeline_test_log.txt" im gleichen Ordner
log_file = "pipeline_test_log.txt"

logging.basicConfig(
    level=logging.DEBUG, # Zeichnet ALLES auf (Infos, Warnungen, Fehler, Debug-Werte)
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    handlers=[
        logging.FileHandler(log_file, mode='w'), # mode='w' überschreibt das Log bei jedem Neustart
        logging.StreamHandler()                  # Gibt es parallel in der Konsole aus
    ]
)

logging.info("=== START PIPELINE TEST RUN ===")


#****INPUT JSON FILE PATH****#
INPUT_JSON = "input.json" 
#****************************#

##################################################
# Parameters are loaded from the input json file #
##################################################

#Entrez request parameters
REQUEST_LIMIT = 5
SLEEP_TIME = .5
EMAIL = ""
E_API = ""

#Blast parameters
blast_type = ''
local_db_path = ''
tax_include = []
tax_exclude = []
database = 'ref_prok_rep_genomes'
e_val = 10-10
coverage_min = .8
max_hits = 500
max_blast_attempts = 3
blast_search_mult_factor = 2
annotate = True
extensive_search = True
reverse_blast = True

#Hit feature detection parameters
margin_limit = 15
max_feature_detect_attempts = 10
feature_search_mult_factor = 4

#Operon assembly parameters
use_reference_threshold = True
ref_threshold_margin = .5
feature_limit = 3
intergenic_limit = 1500

#Other paramters
thread_limit = 10
species_percent_id_limit = 0.45
color_code = {}

# Input / reference configuration
input_records = []
input_record_type = "protein_accession"

reference_database_mode = "remote"   # "remote" or "local"
reference_data = None                # remote: nuccore accession, local: reference file path
reference_genome_accession = ''
reference_assembly_accession = ''
reference_genome_name = ''
reference_cds_select_by = None
reference_cds_select_values = None

#Output parameters
cache_dir = './cache/'
output_dir = './output/{run_id}/'
run_id = str(datetime.date.today()) + '_' + str(uuid.uuid4())

##################################################
##################################################

#Data related to the reference operon
ref_genome_frag = None
ref_features = []

#Reverse BLAST db directories
reverse_blast_root = './reverse_blast/{ref_assembly_accession}'
reference_total_protein = './reverse_blast/{ref_assembly_accession}/{ref_assembly_accession}.fasta'
reference_blast_db = './reverse_blast/{ref_assembly_accession}/{ref_assembly_accession}_blastdb'

#The list of species
species = []

def load_reference_from_local_genbank(genbank_path, cds_select):
    """
    Loads a local nucleotide GenBank file and extracts selected CDS features
    to be used as reference operon genes and BLAST queries.

    Parameters
    ----------
    genbank_path : str
        Path to local GenBank file.
    cds_select : dict
        {
            "by": "locus_tag",
            "values": ["AJMIMBIC_04697", ...]
        }

    Returns
    -------
    reference_record : SeqRecord
        Full GenBank record.
    selected_features : list
        List of selected CDS features.
    query_proteins : list
        List of SeqRecord protein sequences for BLAST.
    """

    if not os.path.exists(genbank_path):
        raise FileNotFoundError(f"GenBank file not found: {genbank_path}")

    reference_record = SeqIO.read(genbank_path, "genbank")

    select_by = cds_select["by"]
    select_values = set(cds_select["values"])

    selected_features = []
    query_proteins = []

    for feature in reference_record.features:

        if feature.type != "CDS":
            continue

        identifier = None

        # Determine which qualifier to use
        if select_by == "locus_tag":
            identifier = feature.qualifiers.get("locus_tag", [None])[0]

        elif select_by == "gene":
            identifier = feature.qualifiers.get("gene", [None])[0]

        elif select_by == "protein_id":
            identifier = feature.qualifiers.get("protein_id", [None])[0]

        else:
            raise ValueError(f"Unsupported cds_select method: {select_by}")

        if identifier in select_values:

            selected_features.append(feature)

            # Extract protein sequence
            translation = feature.qualifiers.get("translation")

            if translation:
                protein_seq = translation[0].replace("\n", "").replace(" ", "")
            else:
                # fallback: translate from nucleotide
                nuc_seq = feature.extract(reference_record.seq)
                protein_seq = str(nuc_seq.translate(to_stop=True))

            # Build SeqRecord for BLAST query
            protein_record = SeqRecord(
                Seq(protein_seq),
                id=identifier,
                description=f"local_ref|{identifier}"
            )

            query_proteins.append(protein_record)

    if len(selected_features) == 0:
        raise ValueError("No matching CDS found in GenBank for provided cds_select values.")

    if len(selected_features) != len(select_values):
        print("WARNING: Not all requested CDS were found in GenBank.")
        print("Requested:", select_values)
        print("Found:", [f.qualifiers.get(select_by, ["-"])[0] for f in selected_features])

    print(f"Loaded {len(selected_features)} CDS from local GenBank reference.")

    return reference_record, selected_features, query_proteins

def search_blast(
    input_records,
    db='ref_prok_rep_genomes',
    max_attempts=3,
    search_mult_factor=2,
    max_hits=50,
    e_cutoff=10E-10,
    tax_incl=[],
    tax_excl=[],
    annotate=True,
    min_cover=None,
    extensive_search=True
):
    '''
    Performs blast search for a set of records. 

    Parameters
    ----------
    input_records : list[string]
        The list of accession numbers (protein IDs OR local locus_tags) to conduct the BLAST search.
    db: string, optional
        The database to use for the BLAST search.
    max_attempts: int, optional
        The BLAST search will be repeated this many times until all hits within the e-value cutoff are returned.
    search_mult_factor: int, optional
        The max_hits will be multiplied by this factor following every incomplete BLAST search.
    max_hits : int, optional
        Starting max number of hits to return.
    e_cutoff : float, optional
        The threshold for the E-value of the hits.
    tax_incl : list[int], optional
        The taxa to include in the BLAST search.
    tax_excl : list[int], optional
        The taxa to exclude in the BLAST search.
    annotate:
        If True return AnnotatedHit objects, otherwise tuple formatted data.
    min_cover: float, optional
        The minimum coverage of the hits.
    extensive_search: bool
        If False, only iterate once.

    Returns
    -------
    hits : list[(input_record, hit_record ,alignment_object)] (if annotate=False)
    annotated_hits: list[AnnotatedHit.object] (if annotate=True)
    '''

    # Holding the max attempts in a temp variable so that it can reset for each input record passed in
    temp_max_hits = max_hits

    # If extensive search is off, the max_attempts is set to 1
    if not extensive_search:
        max_attempts = 1

    print("|~> BLAST search: " + str(input_records) + "...")

    # Check the length of the input_records.
    if len(input_records) < 1:
        raise Exception("Need at least one record to conduct BLAST search.")

    # The final list of BLAST hits.
    return_hits = []

    # Helper: build entrez taxon query string
    def _build_taxon_query(tax_incl_list, tax_excl_list):
        query = ""
        
        # --- INCLUDES ---
        if len(tax_incl_list) == 1:
            query += f"txid{tax_incl_list[0]}[orgn]"
        elif len(tax_incl_list) > 1:
            incl_str = " OR ".join([f"txid{tax}[orgn]" for tax in tax_incl_list])
            query += f"({incl_str})"
            
        # --- EXCLUDES ---
        if len(tax_excl_list) == 1:
            if query != "":
                query += " NOT "
            query += f"txid{tax_excl_list[0]}[orgn]"
        elif len(tax_excl_list) > 1:
            excl_str = " OR ".join([f"txid{tax}[orgn]" for tax in tax_excl_list])
            if query != "":
                query += " NOT "
            query += f"({excl_str})"
            
        return query

    # Helper: get query sequence either from local reference (locus_tag) or from NCBI protein efetch (protein accession)
    def _get_query_seqrecord(input_record):
        """
        Returns a SeqRecord with the amino-acid sequence for the query.
        - If input is raw sequence: returns it directly.
        - If local genbank reference mode: uses ref_features to find aa_sequence
        - Otherwise: uses Entrez efetch on db=protein
        """
        # Check if user explicitly passed a sequence based on the input_record_type
        if input_record_type == "translation":
            print("DEBUG: Input appears to be a raw amino acid sequence.")
            # Generate a dummy ID
            dummy_id = "SEQ_" + str(input_record)[:10] 
            return SeqRecord(Seq(input_record), id=dummy_id, description="user_provided_translation")
        
        # Detect local reference mode by presence of reference_data AND ref_features
        local_mode = (
            ('reference_database_mode' in globals())
            and (reference_database_mode == "local")
            and ('ref_features' in globals())
            and (len(ref_features) > 0)
        )

        if local_mode:
            # input_record is expected to be a locus_tag (as you set in JSON),
            # but we allow matching by locus_tag OR protein_accession just in case.
            feat = None
            for f in ref_features:
                if getattr(f, "locus_tag", None) == input_record or getattr(f, "protein_accession", None) == input_record:
                    feat = f
                    break

            if feat is None:
                raise ValueError(f"Local reference mode: could not find feature for '{input_record}' in ref_features.")

            aa = getattr(feat, "aa_sequence", None)
            if not aa:
                raise ValueError(f"Local reference mode: feature '{input_record}' has no aa_sequence (translation missing).")

            print("DEBUG QUERY SEQ LENGTH:", len(aa))

            return SeqRecord(Seq(aa), id=str(input_record), description="local_genbank_query")

        # --- OLD behavior: input_record must be a protein accession resolvable via Entrez ---
        print("\t|~> Getting protein record for " + str(input_record) + "...")
        handle = None
        for i in range(REQUEST_LIMIT):
            try:
                handle = Entrez.efetch("protein", id=input_record, rettype="fasta", retmode="text")
                time.sleep(SLEEP_TIME)
                break
            except Exception:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                time.sleep(SLEEP_TIME)
                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

        if handle is None:
            raise ValueError(f"Could not fetch protein FASTA for '{input_record}' from NCBI after {REQUEST_LIMIT} attempts.")

        print("\t|~> Getting protein sequence for " + str(input_record) + "...")
        try:
            return SeqIO.read(handle, "fasta")
        except Exception as e:
            raise ValueError(f"Failed to parse FASTA for '{input_record}'. NCBI returned empty/invalid FASTA. ({e})")

    # Gets the accession numbers for all the hits in the BLAST search and appends return_hits with every unique record.
    for input_record in input_records:

        # Reset the max_hits for each input record
        max_hits = temp_max_hits

        # Fetch query sequence (local or Entrez)
        input_seq = _get_query_seqrecord(input_record)

        # Keeps track of attempts
        current_attempts = 0
        blast_complete = False

        # Continue repeating BLAST searches until:
        # - max_attempts is reached
        # - or all hits within the e-value cutoff are returned.
        while not blast_complete:

            if current_attempts == max_attempts:
                print("\t\t|~> Max number of BLAST search attempts has been reached.")
                blast_complete = True
                continue

            print("\t|~> Performing BLAST search " + str(current_attempts + 1))

            blast_records = None

            # Build entrez query if needed
            use_taxon = (len(tax_incl) > 0 or len(tax_excl) > 0)
            taxon_query = _build_taxon_query(tax_incl, tax_excl) if use_taxon else None

            # Send BLAST request
            for i in range(REQUEST_LIMIT):
                try:
                    if use_taxon:
                        result_handle = NCBIWWW.qblast(
                            "tblastn",
                            db,
                            input_seq.format('fasta'),
                            entrez_query=taxon_query,
                            expect=e_cutoff,
                            hitlist_size=max_hits
                        )
                    else:
                        result_handle = NCBIWWW.qblast(
                            "tblastn",
                            db,
                            input_seq.format('fasta'),
                            expect=e_cutoff,
                            hitlist_size=max_hits
                        )

                    print("\t\t|~> Getting BLAST result records")
                    blast_records = list(NCBIXML.parse(result_handle))
                    time.sleep(SLEEP_TIME)
                    break

                except Exception:
                    print("\t\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    time.sleep(SLEEP_TIME)
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

            # If BLAST failed entirely
            if blast_records is None or len(blast_records) == 0:
                print("\t|~> BLAST search returned no parsable records. Reattempting...")
                current_attempts += 1
                continue

            # If number of returned hits is zero, continue
            if len(blast_records[0].alignments) == 0:
                print("\t|~> BLAST search returned no hits. Reattempting...")
                current_attempts += 1
                continue

            # If fewer hits than requested, likely complete
            if len(blast_records[0].alignments) < max_hits:
                blast_complete = True
                print("\t|~> BLAST search was successful")
                continue

            # Check E-value of last hit
            if blast_records[0].alignments[-1].hsps[0].expect > e_cutoff:
                blast_complete = True
                print("\t|~> BLAST search was successful")
                continue
            else:
                print("\t|~> BLAST search was not complete. Increasing max_hits...")
                max_hits = max_hits * search_mult_factor
                print("\t|~> Max hits: " + str(max_hits) + ". Reattempting...")
                current_attempts += 1

        # Extract hits
        print("\t|~> Extracting hits from BLAST results...")
        for record in blast_records[0].alignments:

            current_hit_def = re.sub('[^A-Za-z0-9]+', '_', record.hit_def)
            curr_hit_rec = record.hit_id.split('|')[-2]
            print("\t\t|~> Analyzing hit " + str(curr_hit_rec))

            for hit in record.hsps:

                # If long sequence is passed by the user we have to create a short query ID
                query_id_for_hit = str(input_record)
                if len(query_id_for_hit) > 20:
                    query_id_for_hit = "SEQ_" + query_id_for_hit[:10]

                # Initiate AnnotatedHit if requested
                if annotate:
                    a_hit = AnnotatedHit(
                        query_accession=query_id_for_hit,
                        hit_accession=curr_hit_rec,
                        genome_fragment_name=current_hit_def,
                        align_start=hit.sbjct_start,
                        alignment_seq=hit.sbjct,
                        align_end=hit.sbjct_end,
                        strand=hit.frame[1],
                        percent_identity=(hit.identities / hit.align_length),
                        req_limit=REQUEST_LIMIT,
                        sleep_time=SLEEP_TIME
                    )

                # Coverage filter
                if min_cover:
                    cov = (hit.query_end - hit.query_start + 1) / (len(input_seq.seq))
                    print('\t\t\tCoverage value: ' + str(cov))

                    if cov >= min_cover:
                        if annotate:
                            return_hits.append(a_hit)
                        else:
                            return_hits.append((query_id_for_hit, curr_hit_rec, record))
                    else:
                        print("\t\t|~> Hit did not meet coverage requirement: " + str(curr_hit_rec))
                        print('\t\t\tCoverage value: ' + str(cov))

                else:
                    # Append unique hits
                    if annotate:
                        if len(return_hits) == 0:
                            print("\t\t|~> Adding first hit: " + str(curr_hit_rec))
                            return_hits.append(a_hit)
                        elif not (a_hit in return_hits):
                            print("\t\t|~> Adding hit: " + str(curr_hit_rec))
                            return_hits.append(a_hit)
                    else:
                        if len(return_hits) == 0:
                            print("\t\t|~> Adding first hit: " + str(curr_hit_rec))
                            return_hits.append((query_id_for_hit, curr_hit_rec, record))
                        elif not (curr_hit_rec in list(zip(*return_hits))[1]):
                            print("\t\t|~> Adding hit: " + str(curr_hit_rec))
                            return_hits.append((query_id_for_hit, curr_hit_rec, record))

    print("\t|~> Returning " + str(len(return_hits)) + " unique hits")
    return return_hits

def local_blast_search(input_record, db_path, e_cutoff=10-10, min_cover=None):
    '''
    Completes a local blast search instead of conducting remotely.

    Parameters
    ----------
    input_record: str
        The accession for the query protein record
    db_path: str
        The file path to the local BLAST database
    e_cutoff: float
        The maximum e-value to cutoff the blast hits.
    min_cover: float
        The minimum coverage for each accepted hit
    
    Returns
    -------
    annotated_hits: list[AnnotateHit.object]
        A list of AnnotatedHit.objects that hold metadata for each of the BLAST hits.
    '''

    print('Local BLAST: ' + str(input_record))

    # Validate DB prefix early (Biopython will otherwise run `tblastn -db` with an empty value and only print USAGE)
    if db_path is None or str(db_path).strip() == '':
        raise ValueError(
            "blast_type is 'local' but 'blast.local_db_path' is empty in your input JSON. "
            "Set it to the BLAST database prefix (path WITHOUT extension), e.g. './db/my_nucl_db'."
        )

    # `-db` takes a prefix, not a single file; check that some DB files exist for that prefix
    db_prefix = str(db_path).strip()
    if len(glob.glob(db_prefix + '.*')) == 0:
        raise FileNotFoundError(
            f"Local BLAST database prefix not found: '{db_prefix}'. "
            f"Expected files like '{db_prefix}.nin'/'{db_prefix}.nsq' (nucl) or '{db_prefix}.pin'/'{db_prefix}.psq' (prot)."
        )

    #Get the fasta record for the input record
    global input_record_type
    fasta_record = None

    if input_record_type == "translation":
        print("\t|~> Input is a raw amino acid sequence. Skipping NCBI download.")
        # Wir bauen uns das FASTA-Format einfach selbst
        dummy_id = "SEQ_" + str(input_record)[:10]
        fasta_record = f">{dummy_id}\n{input_record}\n"
    else:
        # Get the fasta record for the input record via NCBI
        for i in range(REQUEST_LIMIT):
            try:
                handle = Entrez.efetch(db='protein', id=input_record, retmode='fasta', rettype='fasta')        
                fasta_record = handle.read()
                time.sleep( SLEEP_TIME)
                break
            except Exception as e:
                print(f"\t\tNCBI exception raised on attempt {i+1} for {input_record}: {e}\n\t\treattempting now...")
                if i == (REQUEST_LIMIT - 1):
                    print(f"\t\tCould not download record after {REQUEST_LIMIT} attempts")

    #Check if the fasta record was pulled successfully
    if fasta_record == None:
        print('\t\tFasta record could not be downloaded for ' + str(input_record))
        return None
    
    os.makedirs('./local_blast_bin/', exist_ok=True)
    
    #Write the FASTA record to a temporary input file for the local blast search
    record_fasta_file = './local_blast_bin/temp_in.fasta'
    with open(record_fasta_file, 'w') as file:
        file.write(fasta_record)
    
    #Get the query length that will be used to calculate the coverage later
    record_fasta = SeqIO.read(open(record_fasta_file,'r'),'fasta')
    query_length = len(record_fasta.seq)


    print('Downloaded query sequence: ' + str(input_record))

    #Conduct the local BLAST search
    blast_command = NcbitblastnCommandline(query=record_fasta_file, db=db_path, evalue=e_cutoff, outfmt=5, out="./local_blast_bin/out.xml")
    
    try:
        blast_command()
    except Exception as e:
        raise RuntimeError(f"Local BLAST command failed! Make sure NCBI BLAST+ is installed on your computer. Error: {e}")
    
    # Parse the BLAST results
    with open('./local_blast_bin/out.xml', 'r') as out_handle:
        blast_records = list(NCBIXML.parse(out_handle))

    print('Parsing through local BLAST results ' + str(input_record) + '...')

    #List of annotated hits to return
    return_hits = []

    print("\t|~> Extracting hits from BLAST results...")

    if not blast_records or len(blast_records[0].alignments) == 0:
        print("\t\t|~> No hits found in local database.")
        return return_hits

    for record in blast_records[0].alignments:
        hit_def_parts = record.hit_def.split(' ')
        current_hit_def = re.sub('[^A-Za-z0-9]+', '_', hit_def_parts[1] if len(hit_def_parts) > 1 else hit_def_parts[0])
        curr_hit_rec = hit_def_parts[0]
        
        print("\t\t|~> Analyzing hit " + str(curr_hit_rec))
        
        #Iterate through the hits
        for hit in record.hsps:

            #Initiates a AnnotatedHit object if set by the parameters.
            a_hit = AnnotatedHit(
                query_accession=input_record, 
                hit_accession=curr_hit_rec, 
                genome_fragment_name=current_hit_def, 
                align_start=hit.sbjct_start, 
                alignment_seq=hit.sbjct, 
                align_end=hit.sbjct_end, 
                strand=hit.frame[1], 
                percent_identity=(hit.identities/hit.align_length), 
                req_limit=REQUEST_LIMIT, 
                sleep_time=SLEEP_TIME
            )

            if min_cover == None:
                if annotate:
                    return_hits.append(a_hit)
                else:
                    return_hits.append((input_record, curr_hit_rec, record))
                continue

            #Calculate the coverage for the current hit                  
            cov = (hit.query_end - hit.query_start + 1) / (query_length)
            print('\t\t\tCoverage value: ' + str(cov))
            
            if cov >= min_cover:
                if annotate:
                    return_hits.append(a_hit)
                else:
                    return_hits.append((input_record, curr_hit_rec ,record))
            else:
                print("\t\t|~> Hit did not meet coverage requirement: " + str(curr_hit_rec))
                
    print("\t|~> Returning " + str(len(return_hits)) + " unique hits")
    return return_hits

def load_input_file(filename):
    '''
    Loads all the paramters from the input JSON.

    Parameters
    ----------
    filename: string
        The input file located inside the /input directory that is in the same directory as this .py script.
    '''

    # Open the JSON file
    file_reader = json.load(open(filename))

    # ----------------------------
    # Entrez request parameters
    # ----------------------------
    global REQUEST_LIMIT
    REQUEST_LIMIT = file_reader['entrez'][0]['request_limit']

    global SLEEP_TIME
    SLEEP_TIME = file_reader['entrez'][0]['sleep_time']

    global EMAIL
    EMAIL = file_reader['entrez'][0]['email']

    global E_API
    E_API = file_reader['entrez'][0]['api_key']

    # ----------------------------
    # Blast parameters (with safe defaults for backwards compatibility)
    # ----------------------------
    global blast_type
    blast_type = file_reader['blast'][0].get('blast_type', 'remote')

    global local_db_path
    local_db_path = file_reader['blast'][0].get('local_db_path', '')

    global tax_include
    tax_include = file_reader['blast'][0].get('tax_include', [])

    global tax_exclude
    tax_exclude = file_reader['blast'][0].get('tax_exclude', [])

    global database
    database = file_reader['blast'][0].get('database', 'ref_prok_rep_genomes')

    global e_val
    e_val = file_reader['blast'][0].get('e_val', 1e-10)

    global coverage_min
    coverage_min = file_reader['blast'][0].get('coverage_min', 0.5)

    global max_hits
    max_hits = file_reader['blast'][0].get('max_hits', 20000)

    global max_blast_attempts
    max_blast_attempts = file_reader['blast'][0].get('max_attempts', 3)

    global blast_search_mult_factor
    blast_search_mult_factor = file_reader['blast'][0].get('search_mult_factor', 2)

    global annotate
    annotate = file_reader['blast'][0].get('annotate', True)

    global extensive_search
    extensive_search = file_reader['blast'][0].get('extensive_search', True)

    global reverse_blast
    reverse_blast = file_reader['blast'][0].get('reverse_blast', False)

    # ----------------------------
    # Hit feature detection parameters
    # ----------------------------
    global margin_limit
    margin_limit = file_reader['hit_feature_detection'][0]['margin_limit']

    global max_feature_detect_attempts
    max_feature_detect_attempts = file_reader['hit_feature_detection'][0]['max_attempts']

    global feature_search_mult_factor
    feature_search_mult_factor = file_reader['hit_feature_detection'][0]['mult_factor']

    # ----------------------------
    # Operon assembly parameters
    # ----------------------------
    global feature_limit
    feature_limit = file_reader['operon_assembly'][0]['feature_limit']

    global intergenic_limit
    intergenic_limit = file_reader['operon_assembly'][0]['intergenic_limit']

    global use_reference_threshold
    use_reference_threshold = file_reader['operon_assembly'][0]['use_ref_limit']

    global ref_threshold_margin
    ref_threshold_margin = file_reader['operon_assembly'][0]['ref_limit_margin']

    # ----------------------------
    # Other parameters
    # ----------------------------
    global thread_limit
    thread_limit = file_reader.get('thread_limit', 1)

    global species_percent_id_limit
    species_percent_id_limit = file_reader.get('species_percent_id_limit', 0.0)

    global color_code
    # allow missing color_code for quick tests
    color_code = file_reader.get('color_code', [{}])[0]

    # ----------------------------
    # Reference / Input records
    # ----------------------------
        # ----------------------------
    # Reference / Input records
    # ----------------------------
    global input_records
    global input_record_type
    global reference_database_mode
    global reference_data
    global reference_genome_accession
    global reference_assembly_accession
    global reference_genome_name
    global reference_cds_select_by
    global reference_cds_select_values

    # defaults
    input_records = []
    input_record_type = "protein_accession"
    reference_database_mode = "remote"
    reference_data = None
    reference_genome_accession = ''
    reference_assembly_accession = ''
    reference_genome_name = ''
    reference_cds_select_by = None
    reference_cds_select_values = None

    # New unified schema
    if 'reference' in file_reader:
        ref_cfg = file_reader['reference']

        reference_database_mode = ref_cfg.get('database_mode', None)
        reference_data = ref_cfg.get('data', None)

        if reference_database_mode not in ("local", "remote"):
            raise KeyError("reference.database_mode must be either 'local' or 'remote'")

        if not reference_data:
            raise KeyError("reference.data is missing or empty")

        # Input records are now unified as:
        # "input_records": {"type": "...", "values": [...]}
        if 'input_records' not in file_reader:
            raise KeyError("Missing 'input_records' block")

        input_cfg = file_reader['input_records']
        input_record_type = input_cfg.get('type', None)
        input_records = input_cfg.get('values', None)

        if not input_record_type:
            raise KeyError("input_records.type is missing")
        if not input_records:
            raise KeyError("input_records.values is missing or empty")

        # REMOTE MODE
        if reference_database_mode == "remote":
            reference_genome_accession = reference_data
            reference_genome_name = ref_cfg.get('reference_genome_name', None)

            if not reference_genome_name:
                raise KeyError("reference.reference_genome_name is required in remote mode")

            # optional; can be derived later from nuccore accession
            reference_assembly_accession = ref_cfg.get('reference_genome_assembly', '')

        # LOCAL MODE
        elif reference_database_mode == "local":
            # reference_data is the path to the anchor local GenBank file
            # input_records are expected to be locus_tags / protein_ids depending on input_record_type
            reference_cds_select_by = input_record_type
            reference_cds_select_values = list(input_records)

            # keep remote-only values empty
            reference_genome_accession = ''
            reference_assembly_accession = ''
            reference_genome_name = ''

    else:
        # ----------------------------
        # Old schema fallback (optional; remove later if desired)
        # ----------------------------
        if 'input_records' not in file_reader:
            raise KeyError("Missing 'input_records' and no 'reference' block found")

        input_records = file_reader['input_records']
        input_record_type = "protein_accession"

        reference_database_mode = "remote"
        reference_data = file_reader['reference_genome_accession']
        reference_genome_accession = file_reader['reference_genome_accession']
        reference_assembly_accession = file_reader.get('reference_genome_assembly', '')
        reference_genome_name = file_reader['reference_genome_name']

    if len(input_records) < 2:
        raise ValueError("Invalid input: Operon conservation analysis requires at least two genes to calculate structural pairs. Only one was provided.")
    
    # ----------------------------
    # Output parameters
    # ----------------------------
    global cache_dir
    cache_dir = file_reader.get('cache_dir', './cache/')
    if cache_dir is None or cache_dir == "":
        cache_dir = './cache/'
    cache_dir = str(cache_dir).replace('\\', '/')
    if not cache_dir.endswith('/'):
        cache_dir = cache_dir + '/'

    global output_dir
    output_dir = file_reader.get('output_dir', 'output/')
    if output_dir is None or output_dir == "":
        output_dir = 'output/'
    output_dir = str(output_dir).replace('\\', '/')
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    
    os.makedirs(cache_dir, exist_ok=True)
    

def write_all_out(species_list, query_accessions, output_path):
    '''
    Writes the results into two separate CSV files: 
    1. Summary: Best operon and genomic totals per species.
    2. Detailed: Every single operon fragment found.
    '''
    import csv
    import os

    summary_file = output_path.replace(".csv", "_summary.csv")

    # WRITE SUMMARY FILE
    with open(summary_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        header = [
            'Species', 
            'TaxID', 
            'Assembly_Accession', 
            'Best_Operon_Nucleotide_ID', 
            'Best_Operon_SIM', 
            'Total_Genomic_SIM', 
            'Best_Operon_Avg_AAI', 
            'Global_Total_Avg_AAI'
        ]
        
        for q_acc in query_accessions:
            header.append(f'Best_Op_{q_acc}_AAI')
            
        for q_acc in query_accessions:
            header.append(f'Total_Count_{q_acc}')
            
        writer.writerow(header)

        for sp in species_list:
            row = []
            row.append(sp.species_name)
            row.append(getattr(sp, 'taxid', 'N/A'))
            row.append(getattr(sp, 'assembly_accession', 'N/A'))
            
            if sp.best_operon:
                row.append(sp.best_operon.genome_accession)
            else:
                row.append('N/A')

            if sp.best_operon:
                row.append(f"{sp.best_operon.local_sim * 100:.1f}%")
                row.append(f"{sp.sim_score * 100:.1f}%")
                row.append(f"{sp.best_operon.local_aai * 100:.1f}%")
            else:
                row.extend(["0.0%", f"{sp.sim_score * 100:.1f}%", "0.0%"])

            global_aai = getattr(sp, 'total_avg_aai_calculated', 0.0)
            row.append(f"{global_aai * 100:.2f}%")

            for q_acc in query_accessions:
                aai_val = "N/A"
                if sp.best_operon:
                    for feat in sp.best_operon.features:
                        if isinstance(feat, AnnotatedHit) and feat.query_accession == q_acc:
                            aai_val = f"{feat.percent_identity * 100:.1f}%"
                            break
                row.append(aai_val)

            for q_acc in query_accessions:
                count = sp.query_hits_counts.get(q_acc, 0)
                row.append(count)

            writer.writerow(row)
        
def append_detailed_out(sp, query_accessions, detailed_path):
    '''
    Appends the operon fragments of a single species to the detailed CSV.
    '''
    import csv
    import os
    from features import AnnotatedHit # Sicherstellen, dass der Import vorhanden ist

    file_exists = os.path.isfile(detailed_path)
    
    with open(detailed_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        if not file_exists:
            header = ['Species', 'Nucleotide Accession', 'Fragment_SIM', 'Fragment_Avg_AAI']
            for q_acc in query_accessions:
                header.extend([
                    f'{q_acc}_AAI', 
                    f'{q_acc}_Accession', 
                    f'{q_acc}_Start', 
                    f'{q_acc}_Stop', 
                    f'{q_acc}_Strand'
                ])
            writer.writerow(header)

        for frag in sp.genome_fragments:
            for operon in frag.operons:
                row = [
                    sp.species_name, 
                    operon.genome_accession, 
                    f"{operon.local_sim * 100:.1f}%", 
                    f"{operon.local_aai * 100:.1f}%"
                ]
                
                for q_acc in query_accessions:
                    found_hit = None
                    for feat in operon.features:
                        if isinstance(feat, AnnotatedHit) and feat.query_accession == q_acc:
                            found_hit = feat
                            break
                    
                    if found_hit:
                        row.extend([
                            f"{found_hit.percent_identity * 100:.1f}%", 
                            found_hit.protein_accession, 
                            found_hit.five_end, 
                            found_hit.three_end, 
                            found_hit.strand
                        ])
                    else:
                        row.extend(['N/A', 'N/A', 'N/A', 'N/A', 'N/A'])
                
                writer.writerow(row)

def run_itol_pipeline():
        """
        Run the iTOL tree/dataset generation script after output.csv has been written.

        The iTOL files are written into a folder named 'itol' next to the main output.csv.
        """

        global output_dir
        global EMAIL
        global E_API

        # Path to the CSV written by write_all_out()
        output_csv = os.path.join(output_dir, "output_summary.csv")

        if not os.path.exists(output_csv):
            print("iTOL pipeline skipped: output_summary.csv was not found.")
            return

        # Create output_dir/itol
        itol_outdir = os.path.join(output_dir, "itol")
        os.makedirs(itol_outdir, exist_ok=True)

        # Path to the iTOL script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        itol_script = os.path.join(script_dir, "itol_pipeline", "build_itol_tree_and_datasets.py")

        if not os.path.exists(itol_script):
            print("iTOL pipeline skipped: build_itol_tree_and_datasets.py was not found.")
            print("Expected at:", itol_script)
            return

        cmd = [
            sys.executable,
            itol_script,
            output_csv,
            "--email", EMAIL,
            "--outdir", itol_outdir
        ]

        # Only pass API key if available
        if E_API and str(E_API).strip() != "":
            cmd.extend(["--api-key", E_API])

        print("Running iTOL pipeline...")
        print("Command:", " ".join(cmd))

        proc = subprocess.run(cmd, capture_output=True, text=True)

        if proc.returncode != 0:
            print("iTOL pipeline failed.")
            print(proc.stdout)
            print(proc.stderr)
        else:
            print("iTOL pipeline finished successfully.")
            print(proc.stdout)

def fetch_nuccore_metadata(curr_acc):
    '''
    Fetches the assembly accession, species name, and taxid for a single nuccore accession.
    Returns a lightweight dictionary to avoid holding heavy objects in memory.
    '''
    metadata = {
        'assembly_accession': curr_acc, # Fallback
        'taxid': '-1',
        'species_name': None
    }
    
    logging.info(f"Downloading metadata records for {curr_acc}...")

    records = None
    for i in range(REQUEST_LIMIT):
        try:
            handle = Entrez.efetch(db="nuccore", id=curr_acc, rettype='gb', retmode='XML')
            records = list(Entrez.read(handle, 'xml'))
            time.sleep(SLEEP_TIME)
            break
        except Exception as e:
            logging.warning(f"NCBI efetch exception on attempt {i+1}/{REQUEST_LIMIT} for {curr_acc}: {e}")
            time.sleep(SLEEP_TIME * 2)

    if not records:
        logging.error(f"Could not fetch data for {curr_acc}. Skipping metadata extraction.")
        return metadata

    record = records[0]
    assembly_found = False

    # 1) Assembly via GBSeq_xrefs
    for info in record.get('GBSeq_xrefs', []):
        if info.get('GBXref_dbname') == 'Assembly':
            metadata['assembly_accession'] = info.get('GBXref_id')
            assembly_found = True
            break

    # 2) Fallback: nuccore -> assembly via elink
    if not assembly_found:
        link_records = None
        for attempt in range(REQUEST_LIMIT):
            try:
                link_handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=curr_acc)
                link_records = Entrez.read(link_handle)
                time.sleep(SLEEP_TIME)
                break
            except Exception:
                time.sleep(SLEEP_TIME)
        
        if link_records:
            asm_id = None
            for lr in link_records:
                for lset in lr.get("LinkSetDb", []):
                    for link in lset.get("Link", []):
                        asm_id = link.get("Id")
                        if asm_id: break
                    if asm_id: break
                if asm_id: break

            if asm_id:
                for attempt in range(REQUEST_LIMIT):
                    try:
                        sum_handle = Entrez.esummary(db="assembly", id=asm_id, retmode="xml")
                        sum_records = Entrez.read(sum_handle)
                        time.sleep(SLEEP_TIME)
                        try:
                            docsum = sum_records['DocumentSummarySet']['DocumentSummary'][0]
                            asm_acc = docsum.get('AssemblyAccession', None)
                            if asm_acc:
                                metadata['assembly_accession'] = asm_acc
                                assembly_found = True
                        except Exception:
                            pass
                        break
                    except Exception:
                        time.sleep(SLEEP_TIME)

    # 3) TaxID + organism name from source feature
    for feature in record.get('GBSeq_feature-table', []):
        if feature.get('GBFeature_key') == 'source':
            for qual in feature.get('GBFeature_quals', []):
                if qual.get('GBQualifier_name') == 'db_xref':
                    val = qual.get('GBQualifier_value', '')
                    if val.startswith('taxon:'):
                        metadata['taxid'] = val.split(':', 1)[1]
                if qual.get('GBQualifier_name') == 'organism':
                    metadata['species_name'] = qual.get('GBQualifier_value')

    return metadata

def process_frag(fragment, lock):
    '''
    Will complete the GenomeFragment object that is passed in by:
        1. Fetching all the features for the fragment
        2. Assign the features for the hit
        3. Assemble the operon
    
    Note: This function was implemented so that the fragments could be processed on multiple threads. 
    
    Parameters
    ----------
    fragment: GenomeFragment object
        The fragment to be processed.
    lock: ThreadLock
        Needed to print out to screen in sync
    '''

    try:
        fragment.fetch_record()
    except Exception as e:
        print(f"[ERROR] fetch_record failed for {fragment.genome_fragment_accession}: {e}")
        return

    fragment.fetch_features()
    fragment.fetch_hit_features(margin_limit=margin_limit, max_attempts=max_feature_detect_attempts, mult_factor=feature_search_mult_factor)
    fragment.purge_hits()
    fragment.assemble_operons(feature_limit=feature_limit, intergenic_limit=intergenic_limit)
    lock.acquire()
    tqdm.write("Completed:\n" + str(fragment) + "-"*50)
    lock.release()

    #Clear up memory by deleting the full features list
    del fragment.full_record
    fragment.full_record = None

def align_input_records_to_biology():
        '''
        Sorts ref_features biologically (5' -> 3') depending on the strand
        and overwrites input_records so the rest of the pipeline uses the true physical order.
        '''
        global ref_features
        global input_records
        global input_record_type

        if not ref_features or len(ref_features) == 0:
            return

        # Determine the strand from the first feature (sanity check ensures all are identical)
        strand = ref_features[0].strand

        if str(strand) == '-1' or str(strand) == '-':
            # MINUS STRAND: Transcription reads backwards (high to low coordinates)
            # We sort descending (reverse=True) using the highest coordinate of each gene
            ref_features = sorted(ref_features, key=lambda x: max(x.five_end, x.three_end), reverse=True)
            logging.info("Sorting reference operon for MINUS strand (descending coordinates).")
        else:
            # PLUS STRAND: Transcription reads forwards (low to high coordinates)
            # We sort ascending using the lowest coordinate of each gene
            ref_features = sorted(ref_features, key=lambda x: min(x.five_end, x.three_end))
            logging.info("Sorting reference operon for PLUS strand (ascending coordinates).")

        # Re-build input_records to match this physical sequence perfectly
        sorted_inputs = []
        
        # We use a set to avoid appending duplicates if a gene was fetched weirdly
        seen = set()
        
        for feat in ref_features:
            target_id = None 
            if input_record_type == "locus_tag" and getattr(feat, 'locus_tag', None) in input_records:
                target_id = getattr(feat, 'protein_accession', None)
                if target_id in [None, 'None', '']:
                    target_id = getattr(feat, 'locus_tag', None)
                
            elif input_record_type == "protein_accession" and getattr(feat, 'protein_accession', None) in input_records:
                target_id = feat.protein_accession
                
            elif input_record_type == "translation" and getattr(feat, 'aa_sequence', None):
                for acc in input_records:
                    if acc in feat.aa_sequence or feat.aa_sequence in acc:
                        target_id = getattr(feat, 'protein_accession', getattr(feat, 'locus_tag', None))
                        break
            else:
                # Fallback
                for acc in input_records:
                    if acc == getattr(feat, 'locus_tag', None) or acc == getattr(feat, 'protein_accession', None):
                        target_id = getattr(feat, 'protein_accession', acc)
                        break
                        
            if target_id and target_id not in seen:
                sorted_inputs.append(target_id)
                seen.add(target_id)

        # Overwrite the global list with the newly found Protein IDs
        input_records = sorted_inputs
        
        if input_record_type in ["translation", "locus_tag"]:
            logging.info(f"Input was '{input_record_type}'. Normalizing to 'protein_accession' for downstream pipeline.")
            input_record_type = "protein_accession"

        logging.info(f"Input records successfully re-ordered and normalized: {input_records}")

def process_reference():
    '''
    Pulls information for the reference operon.

    Supports two modes:
      1) Old mode (NCBI): reference_genome_accession + input_records as protein accessions
      2) New mode (local GenBank nucleotide): reference_data + input_records as locus_tags (or other selector later)

    After this function:
      - ref_genome_frag is a GenomeFragment-like object with .all_features populated
      - ref_features contains the selected features that define the "reference operon"
    '''

    global ref_genome_frag
    global ref_features
    global input_records

    global reference_database_mode
    global reference_data
    global input_record_type

    # Local Mode
    if reference_database_mode == "local":
        from features import GenomeFeature

        print("Processing reference from local GenBank file...")
        print("Reference file:", reference_data)

        record = SeqIO.read(reference_data, "genbank")

        ref_genome_frag = GenomeFragment(
            name=record.description if record.description else record.name,
            genome_fragment_accession=record.id,
            req_limit=REQUEST_LIMIT,
            sleep_time=SLEEP_TIME,
            cache_directory=cache_dir
        )

        ref_genome_frag.full_record = record
        ref_genome_frag.fetch_features()

    elif reference_database_mode == "remote":
        print("Processing reference from remote NCBI GenBank record...")
        print("Reference accession:", reference_genome_accession)

        ref_genome_frag = GenomeFragment(
            name=reference_genome_name,
            genome_fragment_accession=reference_genome_accession,
            req_limit=REQUEST_LIMIT,
            sleep_time=SLEEP_TIME,
            cache_directory=cache_dir
        )

        ref_genome_frag.fetch_record()
        ref_genome_frag.fetch_features()

    else:
        raise ValueError(f"Unsupported reference_database_mode: {reference_database_mode}")

    # Select the reference features based on your input_records (here: locus_tags from cds_select.values)
    wanted = set(input_records)
    reference_features = []

    for gf in ref_genome_frag.all_features:
        if input_record_type == "locus_tag" and gf.locus_tag in wanted:
            reference_features.append(gf)
        elif input_record_type == "protein_accession" and gf.protein_accession in wanted:
            reference_features.append(gf)
        elif input_record_type == "translation" and gf.aa_sequence in wanted:
            reference_features.append(gf)
        elif input_record_type not in ("locus_tag", "protein_accession", "translation"):
            # fallback: accept any
            if gf.locus_tag in wanted or gf.protein_accession in wanted or gf.aa_sequence in wanted:
                reference_features.append(gf)

    if len(reference_features) == 0:
        raise ValueError(
            "No reference CDS matched your selector values. "
            "Check reference[0].cds_select.values vs. the GenBank locus_tag/protein_id fields."
        )
    
    ref_features = reference_features

    # SANITY CHECK: Ensure all reference genes are located on the same strand
    strands = set([gf.strand for gf in reference_features])
    logging.info(f"Reference features sorted. Strands detected: {strands}")
    if len(strands) > 1:
        logging.warning("EDGE CASE DETECTED: Reference genes are on mixed strands!")
        raise ValueError(
            "Warning: The specified reference genes are not located on the same strand! "
            "A biological operon is not possible in this configuration. Please check your JSON input."
        )

    align_input_records_to_biology()

def get_reference_intergenic_distance():
    '''
    Determines the maximum intergenic distance in the reference operon. 

    Parameters
    ----------
    None - already parsed from the input JSON

    Returns
    -------
    intrgnc_dist: int
    The maximum distance between two genes in the reference operon.
    '''

    global ref_genome_frag
    global ref_features
    global ref_threshold_margin

    #Check if the data for the reference operon has been fetched. If not, fetch the information
    if ref_genome_frag == None and len(ref_features) == 0:
        process_reference()
    
    #Determine the maximum intergenic distance
    max_intergenic_distance = 0

    for i in range(len(ref_features)-1):
        #Calculate the intergenic distance between the current feature, i, and the next feature, i+1
        temp_ig_dist = ref_features[i].get_intergenic_distance(ref_features[i+1])
        if  temp_ig_dist > max_intergenic_distance:
            max_intergenic_distance = temp_ig_dist
    
    #Set the intergenic limit and considers the margin
    global intergenic_limit
    intergenic_limit = (1 + ref_threshold_margin) * max_intergenic_distance

def calculate_percent_ids(sp):
    '''
    Determines the percent identity of each of the hits in a Species object by aligning the entire feature sequence with the sequence of the reference hit. This replaces the BLAST percent identity.

    Parameters
    ----------
    None

    Returns
    -------
    None
    '''

    global ref_features

    # A list holding the symbols of the traditional amino acids
    trad_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Iterate through all of the features in every fragment in each species
    for frag in sp.genome_fragments:

        # Metadata logic was moved to write_all_out to keep this function focused on AAI calculation

        for feat in frag.hits:
            
            # Pull the sequence for this feature
            feat_seq = feat.aa_sequence

            # Hold the reference sequence
            ref_seq = None
            query_id = str(feat.query_accession).strip()

            for ref_feat in ref_features:

                ref_ids = set()

                # protein_id
                if getattr(ref_feat, "protein_accession", None):
                    ref_ids.add(ref_feat.protein_accession)
                    ref_ids.add(ref_feat.protein_accession.split('.')[0])

                # locus_tag
                if getattr(ref_feat, "locus_tag", None):
                    ref_ids.add(ref_feat.locus_tag)

                # normalize query
                query_variants = {
                    query_id,
                    query_id.split('.')[0]
                }

                if ref_ids.intersection(query_variants):
                    ref_seq = ref_feat.aa_sequence
                    break
            
            # Purge the reference sequence and the hit feature sequence of any non traditional amino acids
            purged_ref_seq = ''
            purged_feat__seq = ''

            if not ref_seq:
                print(f"[WARN] No reference sequence found for query: {feat.query_accession}. Setting percent identity to 0.")
                ref_seq = ''
            
            if feat_seq == None:
                feat_seq = ''

            for c in ref_seq:
                if c in trad_aa:
                    purged_ref_seq = purged_ref_seq + c

            for c in feat_seq:
                if c in trad_aa:
                    purged_feat__seq = purged_feat__seq + c
            
            ref_seq = purged_ref_seq
            feat_seq = purged_feat__seq

            # Perform the pairwise alignment
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            
            if not ref_seq or not feat_seq or len(ref_seq) == 0 or len(feat_seq) == 0 or str(feat_seq) == 'None':
                logging.info(f"    Gene {query_id} Alignment: Skipped due to empty sequence.")
                feat.percent_identity = 0.0
                continue

            # We only calculate the very first (best) alignment
            alignments = aligner.align(ref_seq, feat_seq)

            try:
                align = next(iter(alignments))
            except StopIteration:
                feat.percent_identity = 0.0
                continue

            align = alignments[0]
            
            formatted_alignment = format(align).split("\n")
            ref_seq_aligned = formatted_alignment[0]
            feat_seq_aligned = formatted_alignment[2]
            matches = 0
            
            align_size = len(ref_seq_aligned)
            size_adj = 0
            
            for x in range(align_size):
                if ref_seq_aligned[x] == "-" or feat_seq_aligned[x] == "-":
                    size_adj += 1
                    continue
                
                if ref_seq_aligned[x] == feat_seq_aligned[x]:
                    matches += 1
                    
            if (align_size - size_adj) > 0:
                average_percent_similar = (matches / (align_size - size_adj))
            else:
                average_percent_similar = 0.0
            
            feat.percent_identity = average_percent_similar

def make_reference_blastdb():
    """
    Builds the local protein BLAST database used for reverse BLAST filtering.

    Supports:
      - Remote reference: derive assembly accession from reference_genome_accession if needed,
        then download fasta_cds_aa from nuccore records linked to the assembly
      - Local GenBank nucleotide reference: extract all CDS translations from the local .gbk
        and build DB from them
    """

    import os
    import glob
    import subprocess

    global reference_assembly_accession
    global reference_genome_accession
    global reverse_blast_root
    global reference_total_protein
    global reference_blast_db
    global reference_database_mode
    global reference_data

    sequences_to_write = []

    # -------------------------------------------------------
    # Determine a stable database ID for output paths
    # -------------------------------------------------------
    if reference_database_mode == "local":
        # Local mode: use local reference file basename as technical DB identifier
        db_id = os.path.splitext(os.path.basename(reference_data))[0]

    elif reference_database_mode == "remote":
        # Remote mode: use assembly accession; derive it from nuccore accession if needed
        if not reference_assembly_accession or str(reference_assembly_accession).strip() == "":
            if not reference_genome_accession or str(reference_genome_accession).strip() == "":
                raise ValueError(
                    "reference_assembly_accession is empty and reference_genome_accession is also missing."
                )
            reference_assembly_accession = get_assembly_from_nuccore_accession(
                nuccore_acc=reference_genome_accession
            )

        db_id = reference_assembly_accession

    else:
        raise ValueError(f"Unsupported reference_database_mode: {reference_database_mode}")

    reference_assembly_accession = db_id  # ensure global variable is set for downstream use
    out_dir = reverse_blast_root.format(ref_assembly_accession=db_id)
    os.makedirs(out_dir, exist_ok=True)

    fasta_path = reference_total_protein.format(ref_assembly_accession=db_id)
    db_prefix = reference_blast_db.format(ref_assembly_accession=db_id)

    # If DB already exists, do nothing (fast re-runs)
    existing = glob.glob(db_prefix + ".*")
    if len(existing) > 0:
        print(f"\tReverse-BLAST DB already exists: {db_prefix}")
        return

    # -------------------------------------------------------
    # LOCAL MODE: build DB from local GenBank nucleotide file
    # -------------------------------------------------------
    if reference_database_mode == "local":
        print(f"\tBuilding reverse-BLAST DB from local GenBank: {reference_data}")

        reference_dir = os.path.dirname(reference_data)
        reference_basename = os.path.basename(reference_data)

        # collect all GenBank-like files in the same directory
        gbk_pattern = ["*.gbk", "*.gb", "*.gbff", "*.genbank"]
        local_files = []

        for pattern in gbk_pattern:
            local_files.extend(glob.glob(os.path.join(reference_dir, pattern)))
        
        # de-duplicate files (in case of overlapping patterns) and sort for reproducibility
        local_files = sorted(list(set(local_files)))

        if len(local_files) == 0:
            raise ValueError(f"No GenBank files found in directory: {reference_dir}")
        
        print(f"\tFound {len(local_files)} local GenBank file(s):")

        for fp in local_files:
            print(f"\t - {os.path.basename(fp)}")

        for gbk_file in local_files:
            print(f"\tProcessing local GenBank file: {gbk_file}")

            try:
                for record in SeqIO.parse(gbk_file, "genbank"):
                    for feat in record.features:
                        if feat.type != "CDS":
                            continue

                        q = feat.qualifiers if feat.qualifiers else {}

                        locus_tag = q.get("locus_tag", [None])[0]
                        protein_id = q.get("protein_id", [None])[0]

                        if protein_id in [None, 'None', '']:
                            if locus_tag not in [None, 'None', '']:
                                protein_id = locus_tag
                            else:
                                raise ValueError(
                                    f"\n[CRITICAL ERROR] No 'protein_id' found for CDS feature (locus_tag: '{locus_tag}') in the local file '{gbk_file}'.\n"
                                    f"Please ensure all CDS features in your local GenBank reference have valid 'protein_id' qualifiers, then try again."
                                )

                        translation = q.get("translation", [None])[0]

                        if translation:
                            aa = translation.replace("\n", "").replace(" ", "")
                        else:
                            # fallback translation from nucleotide
                            try:
                                aa = str(feat.extract(record.seq).translate(to_stop=True))
                            except Exception:
                                continue

                        # In local mode, input_records are usually locus_tags
                        rec_id = locus_tag or protein_id
                        if not rec_id:
                            continue

                        sequences_to_write.append(
                            SeqRecord(Seq(aa), id=str(rec_id), description=f"local_ref|locus_tag={locus_tag}|protein_id={protein_id}|record={record.id}")
                        )
            except Exception as e:
                print(f"\t\tWarning: could not parse {gbk_file}: {e}")

        if len(sequences_to_write) == 0:
            raise ValueError("No CDS translations found in any local GenBank files; cannot build reverse BLAST DB.")
        
        # de-duplicate sequences by FASTA ID
        uniq = {}
        for r in sequences_to_write:
            uniq[r.id] = r
        sequences_to_write = list(uniq.values())

    # -------------------------------------------------------
    # REMOTE MODE: assembly -> nuccore -> fasta_cds_aa
    # -------------------------------------------------------
    elif reference_database_mode == "remote":
        asm_acc = reference_assembly_accession
        print(f"\tBuilding reverse-BLAST DB from NCBI assembly: {asm_acc}")

        # 1) Resolve assembly accession -> assembly UID
        asm_term = f"{asm_acc}[Assembly Accession]"
        search_records = None

        for i in range(REQUEST_LIMIT):
            try:
                handle = Entrez.esearch(
                    db="assembly",
                    term=asm_term,
                    retmode="xml",
                    retmax="20"
                )
                search_records = Entrez.read(handle)
                handle.close()
                time.sleep(SLEEP_TIME)
                break
            except Exception as e:
                print(f"\t\tNCBI exception on attempt {i+1} (assembly esearch): {e}")
                time.sleep(SLEEP_TIME)

        if not search_records or len(search_records.get("IdList", [])) == 0:
            raise RuntimeError(f"Could not find assembly UID for {asm_acc}")

        asm_uid = search_records["IdList"][0]
        print(f"Found assembly UID: {asm_uid}")

        # 2) Resolve assembly UID -> linked nuccore IDs
        def extract_linked_ids(elink_records):
            ids = []
            for r in elink_records:
                for ldb in r.get("LinkSetDb", []):
                    for link in ldb.get("Link", []):
                        _id = link.get("Id")
                        if _id:
                            ids.append(_id)
            return ids

        nuccore_uids = []

        if asm_acc.startswith("GCA_"):
            linknames = ["assembly_nuccore_insdc", "assembly_nuccore_refseq", "assembly_nuccore"]
        elif asm_acc.startswith("GCF_"):
            linknames = ["assembly_nuccore_refseq", "assembly_nuccore_insdc", "assembly_nuccore"]
        else:
            linknames = ["assembly_nuccore_refseq", "assembly_nuccore_insdc", "assembly_nuccore"]

        for ln in linknames:
            elink_records = None

            for i in range(REQUEST_LIMIT):
                try:
                    h = Entrez.elink(
                        dbfrom="assembly",
                        db="nuccore",
                        id=asm_uid,
                        linkname=ln,
                        retmode="xml"
                    )
                    elink_records = Entrez.read(h)
                    h.close()
                    time.sleep(SLEEP_TIME)
                    break
                except Exception as e:
                    print(f"\t\tNCBI exception on attempt {i+1} (elink {ln}): {e}")
                    time.sleep(SLEEP_TIME)

            if elink_records:
                nuccore_uids = extract_linked_ids(elink_records)
                if nuccore_uids:
                    print(f"Using linkname '{ln}' -> {len(nuccore_uids)} nuccore IDs")
                    break

        if not nuccore_uids:
            raise RuntimeError(
                f"Could not find any nuccore UIDs linked to assembly {asm_acc} (UID: {asm_uid})"
            )

        seen = set()
        nuccore_uids = [x for x in nuccore_uids if not (x in seen or seen.add(x))]
        print(f"Total unique nuccore UIDs linked to assembly: {len(nuccore_uids)}")

        # Optional debug output: resolve nuccore UIDs -> accession.version
        try:
            handle = Entrez.esummary(
                db="nuccore",
                id=",".join(nuccore_uids),
                retmode="xml"
            )
            nuccore_summary = Entrez.read(handle)
            handle.close()
            time.sleep(SLEEP_TIME)

            nuccore_accessions = []
            for docsum in nuccore_summary:
                acc = docsum.get("AccessionVersion")
                if acc:
                    nuccore_accessions.append(acc)

            print("Nuccore accessions linked to this assembly:")
            for acc in nuccore_accessions:
                print(acc)
        except Exception as e:
            print(f"\tWarning: could not resolve nuccore accessions for logging: {e}")

        # 3) Fetch CDS protein sequences for all nuccore records
        batch_size = 20

        for start in range(0, len(nuccore_uids), batch_size):
            batch = nuccore_uids[start:start + batch_size]

            for i in range(REQUEST_LIMIT):
                try:
                    handle = Entrez.efetch(
                        db="nuccore",
                        id=",".join(batch),
                        rettype="fasta_cds_aa",
                        retmode="text"
                    )
                    records = list(SeqIO.parse(handle, "fasta"))
                    handle.close()
                    sequences_to_write.extend(records)
                    time.sleep(SLEEP_TIME)
                    break
                except Exception as e:
                    print(
                        f"\t\tNCBI exception on attempt {i+1} "
                        f"(efetch batch {start}-{start+batch_size}): {e}"
                    )
                    time.sleep(SLEEP_TIME)

        if len(sequences_to_write) == 0:
            raise RuntimeError(f"No protein sequences fetched for assembly {asm_acc}")

        # de-duplicate by FASTA ID
        uniq = {}
        for r in sequences_to_write:
            uniq[r.id] = r
        sequences_to_write = list(uniq.values())

    # Write FASTA
    print(f"\tWriting reference proteins to: {fasta_path} (n={len(sequences_to_write)})")
    SeqIO.write(sequences_to_write, fasta_path, "fasta")

    # Build BLAST DB
    print(f"\tRunning makeblastdb for prefix: {db_prefix}")
    cmd = ["makeblastdb", "-in", fasta_path, "-out", db_prefix, "-dbtype", "prot"]
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            "makeblastdb failed:\n" + proc.stdout + "\n" + proc.stderr
        )

    # Verify DB files exist
    created = glob.glob(db_prefix + ".*")
    if len(created) == 0:
        raise RuntimeError("makeblastdb reported success but no DB files were created.")

def _strip_version(acc):
    """
    Strip trailing version from accession-like identifiers, e.g.
    NP_414542.1 -> NP_414542
    """
    if acc is None:
        return None
    acc = str(acc).strip()
    if not acc:
        return None
    return re.sub(r"\.\d+$", "", acc)


def _candidate_id_set(value):
    """
    Build a small normalized identifier set from a single value.
    """
    out = set()
    if value is None:
        return out

    s = str(value).strip()
    if not s:
        return out

    out.add(s)
    out.add(_strip_version(s))
    return {x for x in out if x}


def _get_expected_reverse_ids(query_accession):
    """
    Return all acceptable identifiers for the original reference CDS
    corresponding to query_accession.

    This makes reverse BLAST robust across:
    - local mode (locus_tag queries)
    - remote mode (protein accession queries)
    - exact IDs vs versionless IDs
    """
    global ref_features

    expected = set()
    q = str(query_accession).strip()

    for rf in ref_features:
        rf_locus = getattr(rf, "locus_tag", None)
        rf_prot = getattr(rf, "protein_accession", None)

        # A query can match either field depending on mode
        if q == str(rf_locus) or q == str(rf_prot):
            expected |= _candidate_id_set(rf_locus)
            expected |= _candidate_id_set(rf_prot)

    # fallback: at least include the query itself
    if not expected:
        expected |= _candidate_id_set(q)

    return expected


def _extract_ids_from_reverse_hit(alignment):
    """
    Extract as many useful identifiers as possible from a BLAST top hit.

    Works for:
    - local FASTA IDs
    - NCBI fasta_cds_aa headers
    - headers containing [protein_id=...] and [locus_tag=...]
    - ids like lcl|CP164021.1_prot_XJI05057.1_13
    """
    candidates = set()

    hit_id = getattr(alignment, "hit_id", "")
    hit_def = getattr(alignment, "hit_def", "")

    # raw strings
    for token in [hit_id, hit_def]:
        if token:
            candidates |= _candidate_id_set(token)

    # split hit_id on common separators
    for part in re.split(r"[| \t]", hit_id):
        if part:
            candidates |= _candidate_id_set(part)

    # qualifiers in hit_def, e.g. [protein_id=XJI05057.1], [locus_tag=AB2I46_26220]
    m = re.search(r"\[protein_id=([^\]]+)\]", hit_def)
    if m:
        candidates |= _candidate_id_set(m.group(1))

    # locus_tag extraction (robust)
    m = re.search(r"locus_tag=([^\]| ]+)", hit_def)
    if m:
        val = m.group(1)
        candidates |= _candidate_id_set(val)
        candidates.add(val)

    # FASTA IDs from local / remote patterns:
    # lcl|CP164021.1_prot_XJI05057.1_13
    m = re.search(r"_prot_([^_]+(?:\.\d+)?)_", hit_id)
    if m:
        candidates |= _candidate_id_set(m.group(1))

    # local_ref|AJMIMBIC_04697|recordid
    if hit_def:
        first_token = hit_def.split()[0]
        candidates |= _candidate_id_set(first_token)

    return {x for x in candidates if x}


def check_reverse_blast(query_accession, annotated_hit):
    """
    Reverse BLAST: blastp(hit_sequence) vs reference DB.
    Accept if the TOP hit in the reference DB matches the original query_accession
    (or contains it, depending on header style).
    """

    import os
    import subprocess
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio import Seq
    from Bio.Blast import NCBIXML

    global reference_blast_db
    global reference_assembly_accession

    # paths
    temp_out = './reverse_blast/temp_out.xml'
    temp_in  = './reverse_blast/temp_in.fasta'

    # clean alignment sequence
    purged_alignment_seq = annotated_hit.alignment_seq.replace("-", "")

    if not purged_alignment_seq:
        return False
    
    SeqIO.write(SeqRecord(Seq.Seq(purged_alignment_seq), id='temp'), temp_in, 'fasta')



    db_path = reference_blast_db.format(ref_assembly_accession=reference_assembly_accession)

    # run blastp
    cmd = ["blastp", "-query", temp_in, "-db", db_path, "-outfmt", "5", "-out", temp_out]
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        # show BLAST error clearly
        print("BLASTP failed:", proc.stderr.strip())
        return False

    if not os.path.exists(temp_out) or os.path.getsize(temp_out) == 0:
        return False

    # parse xml
    with open(temp_out, "r") as f:
        try:
            res = list(NCBIXML.parse(f))
        except Exception as e:
            print("Reverse BLAST XML parse failed", e)
            return False

    if not res or len(res[0].alignments) == 0:
        return False

    top = res[0].alignments[0]

    expected_ids = _get_expected_reverse_ids(query_accession)
    observed_ids = _extract_ids_from_reverse_hit(top)

    # Debug output: keep this while testing
    print(f"[RB DEBUG] query={query_accession}")
    print(f"[RB DEBUG] expected_ids={sorted(expected_ids)}")
    print(f"[RB DEBUG] top_hit_id={top.hit_id}")
    print(f"[RB DEBUG] top_hit_def={top.hit_def}")
    print(f"[RB DEBUG] observed_ids={sorted(observed_ids)}")

    if expected_ids.intersection(observed_ids):
        return True
    
    return False


def get_assembly_from_nuccore_accession(nuccore_acc: str, sleep_time: float = 0.34):
    """
    Resolve a nuccore accession (e.g. CP164021.1) to its linked assembly accession.
    Strategy:
      1) Try to read the Assembly xref directly from the nuccore GenBank XML record
      2) If not present, use elink(nuccore -> assembly), then esummary(assembly)
    """
    logging.info(f"Querying nuccore record for: {nuccore_acc}")

    link_records = None
    for attempt in range(REQUEST_LIMIT):
        try:
            link_handle = Entrez.elink(
                dbfrom="nuccore",
                db="assembly",
                id=nuccore_acc,
                retmode="xml"
            )
            link_records = Entrez.read(link_handle)
            link_handle.close()
            time.sleep(SLEEP_TIME)
            break
        except Exception as e:
            logging.warning(f"NCBI elink exception on attempt {attempt+1}/{REQUEST_LIMIT} for {nuccore_acc}: {e}")
            time.sleep(SLEEP_TIME * 2)
            if attempt == REQUEST_LIMIT - 1:
                logging.error(f"Could not perform elink for {nuccore_acc} after {REQUEST_LIMIT} attempts.")
                raise RuntimeError(f"NCBI Entrez.elink failed: {e}")

    assembly_uid = None
    for lr in link_records:
        for linksetdb in lr.get("LinkSetDb", []):
            if "nuccore_assembly" in linksetdb.get("LinkName", ""):
                links = linksetdb.get("Link", [])
                if links:
                    assembly_uid = links[0].get("Id")
                    break
        if assembly_uid:
            break

    if not assembly_uid:
        raise ValueError(f"No linked assembly UID found for {nuccore_acc}")

    print(f"Found assembly UID via elink: {assembly_uid}")

    # Convert assembly UID -> assembly accession
    summary_handle = Entrez.esummary(db="assembly", id=assembly_uid, report="full")
    summary = Entrez.read(summary_handle)
    summary_handle.close()
    time.sleep(SLEEP_TIME)

    synonym = summary["DocumentSummarySet"]["DocumentSummary"][0]["Synonym"]

    # Common field for the assembly accession
    if "_" in reference_genome_accession:
        assembly_acc = synonym["RefSeq"]
    else:
        assembly_acc = synonym["Genbank"]
    if not assembly_acc:
        raise ValueError(f"Could not extract AssemblyAccession from esummary for UID {assembly_uid}")

    print(f"Resolved assembly accession: {assembly_acc}")
    return assembly_acc


def calc_operon_cons():
    '''
    The main function in operon_conserve_detect.py. 
    Pipeline:
        1. Get parameters from the input file
        2. Set up the output directory and copy the input file to it
        3. Conduct BLAST search for each of the individual input records and store them all in a list of hits
        4. Group the hits into their respective GenomeFragment objects
        5. For each GenomeFragment, detect the feature for each hit it contains and assemble the operons for each GenomeFragement
        6. Group the GenomeFragments into Species
        7. Filter the Species based of the species_percent_id_limit
            |-> This will remove any species where no BLAST results are greater than the species_percent_id_limit
        8. Calculate the score for each Species
        9. Output all the operons into a CSV file.
    '''

    global species
    global output_dir
    global run_id
    global use_reference_threshold
    global reverse_blast
    global blast_type
    global local_db_path
    
    #Take note of the start time
    start_time = datetime.datetime.now()

    os.makedirs(cache_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    ##Load parameters
    input_file = ''
    if(len(sys.argv) >= 2):
        input_file = 'input/' + sys.argv[1]
        load_input_file(filename=input_file)
    else:
        input_file = 'input/' + INPUT_JSON
        load_input_file(input_file)
    
    #Set Entrez parameters
    Entrez.email = EMAIL
    Entrez.api_key = E_API

    ##Setup the output directory
    output_dir = output_dir.format(run_id=run_id)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir + 'svg/', exist_ok=True)
    
    #Copy the input file to the ouput directory
    shutil.copy(input_file, output_dir)


    #Process reference 
    tqdm.write('Processing reference...')
    make_reference_blastdb()
    process_reference()
    if use_reference_threshold:
        tqdm.write('Getting reference intergenic distance...')
        get_reference_intergenic_distance()
        tqdm.write(str(intergenic_limit) + '\n')

    ##Conduct BLAST search

    #Holds the hits from all of the BLAST searches
    final_hits = []

    align_input_records_to_biology()

    for record in input_records:

        record_hits = []

        #Conduct a remote blast, if requested
        if blast_type == 'remote':

            record_hits = search_blast(
                input_records=[record],
                db=database,
                max_attempts=max_blast_attempts,
                search_mult_factor=blast_search_mult_factor, 
                min_cover=coverage_min, 
                max_hits=max_hits, 
                tax_incl=tax_include, 
                tax_excl=tax_exclude, 
                e_cutoff=e_val, 
                annotate=annotate, 
                extensive_search=extensive_search
            )

        #Conduct a local blast, if requested
        elif blast_type == 'local':
            record_hits = local_blast_search(input_record=record, db_path=local_db_path, e_cutoff=e_val, min_cover=coverage_min)
        
        else:
            print('Invalid BLAST type...exiting.')
            return

        print(f"DEBUG: Hits for {record}: {len(record_hits)}")

        #Filter via reverse BLAST, if requested 
        if(reverse_blast):
            purged_record_hits = []

            for hit in tqdm(record_hits, desc='Reverse Blast'):
                if check_reverse_blast(record, hit):
                    purged_record_hits.append(hit)

            final_hits.extend(purged_record_hits)
        else:
            final_hits.extend(record_hits)

        print('Running total of hits returned: ' + str(len(final_hits)))
        time.sleep(3)

    ## 1. Map Hits to Nucleotide Accessions (Lightweight Dictionary)
    nuccore_to_hits = {}
    for hit in tqdm(final_hits, desc='Mapping Hits to Nucleotides'):
        if not isinstance(hit, AnnotatedHit):
            print("The annotate option needs to be set to true for grouping to work.")
            return
        if hit.genome_accession not in nuccore_to_hits:
            nuccore_to_hits[hit.genome_accession] = []
        nuccore_to_hits[hit.genome_accession].append(hit)

    ## 2. Pre-process metadata for each unique nucleotide (Lightweight)
    nuccore_metadata = {}
    for acc in tqdm(nuccore_to_hits.keys(), desc='Fetching Metadata'):
        nuccore_metadata[acc] = fetch_nuccore_metadata(acc)

    ## 3. Group nucleotide accessions by Assembly (Species)
    assembly_to_nuccore = {}
    for acc, meta in nuccore_metadata.items():
        asm = meta['assembly_accession']
        if asm not in assembly_to_nuccore:
            assembly_to_nuccore[asm] = []
        assembly_to_nuccore[asm].append(acc)

    ## 4. Iterative Lazy-Loading Processing
    # Load, process, and purge species one by one to save RAM
    global species
    species = [] # List for the final, lightweight summary data

    for asm, acc_list in tqdm(assembly_to_nuccore.items(), desc='Processing Species (Lazy Load)'):
        sp = Species(assembly_accession=asm)
        processing_threads = []

        # A) Load only the fragments for THIS species into RAM
        for acc in acc_list:
            meta = nuccore_metadata[acc]
            frag_name = nuccore_to_hits[acc][0].genome_fragment_name

            frag = GenomeFragment(
                name=frag_name, 
                genome_fragment_accession=acc, 
                req_limit=REQUEST_LIMIT, 
                sleep_time=SLEEP_TIME, 
                cache_directory=cache_dir
            )
            frag.taxid = meta['taxid']
            frag.species_name = meta['species_name']
            frag.assembly_accession = asm

            # Add hits to fragment
            for hit in nuccore_to_hits[acc]:
                frag.add_hit(hit)

            sp.add_genome_fragment(frag)

            # Analyze fragment (loads the .gb file, assembles operons)
            lock = threading.Lock()
            temp_thread = threading.Thread(target=process_frag, args=(frag, lock))
            temp_thread.start()
            processing_threads.append(temp_thread)

        # Wait until all fragments for THIS species are processed
        for t in processing_threads:
            t.join()

        # B) Calculate scores for this species
        tqdm.write(f'Calculating IDs and Scores for {asm}...')
        calculate_percent_ids(sp)
        sp.get_query_percent_ids(input_records)

        above_limit = any(val >= species_percent_id_limit for val in sp.query_percent_ids.values())

        if above_limit:
            sp.extract_names()
            sp.extract_taxids()
            sp.extract_genome_accessions()
            sp.measure_sim(input_records)
            sp.total_avg_aai_calculated = sp.calculate_total_avg_aai()

            # C) Write directly to disk (Lazy Writing)
            detailed_path = os.path.join(output_dir, "output_detailed.csv")
            append_detailed_out(sp, input_records, detailed_path)
            
            sp.draw_figure(color_code=color_code, output_dir=output_dir + 'svg/')

            # D) CLEAN UP RAM
            # sp.clean() deletes the 'sp.genome_fragments' list.
            # Since the global 'genome_frags' list is gone, these 
            # massive objects are now truly destroyed by the Python Garbage Collector.
            sp.clean() 
            
            # The lightweight shell (scores, names) goes into the summary list
            species.append(sp) 
        else:
            tqdm.write(f"Did not meet threshold. Purging {asm} from RAM...")
            sp.clean()

    print("Writing output now...")
    final_output_path = os.path.join(output_dir, "output.csv")
    write_all_out(species, input_records, final_output_path)

    if len(species) == 0:
        print("Skipping iTOL pipeline: no species passed filtering, output.csv contains no data rows.")
    else:
        print("Generating iTOL tree and annotation datasets...")
        try:
            run_itol_pipeline()
        except Exception as e:
            print(f"iTOL pipeline failed: {e}")

    #Get time elapsed
    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time
    print("*"*10 + " Finished in " + str(elapsed_time) + "*"*10)


if __name__ == "__main__":
    calc_operon_cons() 