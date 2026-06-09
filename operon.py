'''
@author: ichaudr

'''

from features import AnnotatedHit, GenomeFeature
import itertools

class Operon:

    def __init__(self, genome_fragment_name, genome_accession, genome_features, strand):
        self.features = []
        self.genome_fragment_name = genome_fragment_name
        self.genome_accession =genome_accession
        self.genome_features = genome_features
        self.strand = strand

        self.local_sim = 0.0
        self.local_aai = 0.0

    def calculate_local_stats(self, reference_sequence, allow_permutations=False):
        '''
        Calculates local structural similarity (SIM) and AAI.
        Includes Discovery Mode to find the best possible permutation fit.
        '''
        hits = [f for f in self.features if isinstance(f, AnnotatedHit)]

        # 1. Calculate local AAI
        if hits:
            self.local_aai = sum(float(h.percent_identity) for h in hits) / len(hits)
        else:
            self.local_aai = 0.0

        if not reference_sequence:
            self.local_sim = 0.0
            self.genetic_order = ""
            return
        
        # Generate the biologically correct sequence
        query_sequence = [h.query_accession for h in hits]
        
        # IMPORTANT FOR CSV: Store the normalized Genetic Order here!
        # Since we corrected the strand above, this is always the true transcriptional direction.
        self.genetic_order = "-".join(query_sequence)

        # Helper function for strict calculation of a specific sequence
        def calc_strict_local_sim(ref_seq):
            ref_pairs_temp = [(ref_seq[i], ref_seq[i+1]) for i in range(len(ref_seq)-1)]
            local_pairs_temp = [(query_sequence[i], query_sequence[i+1]) for i in range(len(query_sequence)-1)]
            found = sum(1 for rp in ref_pairs_temp if rp in local_pairs_temp)
            return found / len(ref_pairs_temp) if ref_pairs_temp else 0.0

        # 2. Standard SIM (Strict Reference Order)
        self.local_sim = calc_strict_local_sim(reference_sequence)

        # 3. Discovery Mode (Iterative test of all permutations)
        self.best_perm_local_sim = self.local_sim
        self.best_perm_order = "-".join(reference_sequence)

        if allow_permutations:
            import itertools
            # Test all full permutations (e.g., ABC, ACB, BAC...)
            for perm in itertools.permutations(reference_sequence):
                p_sim = calc_strict_local_sim(perm)
                # If we find a permutation that fits the fragment better, we overwrite it
                if p_sim > self.best_perm_local_sim:
                    self.best_perm_local_sim = p_sim
                    self.best_perm_order = "-".join(perm)
    
    def add_feature(self, feature):
        '''
        Appends feature to the list of features associated with this operon.
        Normalizes the list immediately: 
          - Plus strand (+): sorted ascending (left to right)
          - Minus strand (-): sorted descending (right to left, true biological direction)
        '''
        if feature in self.genome_features:
            self.features.append(feature)
            
            if self.strand == '+':
                self.features = sorted(self.features, key=lambda f: f.five_end)
            else:
                # On the minus strand, transcription goes from high to low coordinates
                self.features = sorted(self.features, key=lambda f: f.five_end, reverse=True)
    
    
    
    def __str__(self):
        try:
            to_return = "OPERON:" + str(self.genome_fragment_name) + "(" + str(self.genome_accession) + ")\nStrand: " + str(self.strand) + "\n"

            for f in self.features:
                if isinstance(f, AnnotatedHit):
                    to_return = to_return + "\tOriginal Query: " + f.query_accession + "\tHit Accession: " + f.protein_accession + "\tLocus Tag: " + f.locus_tag + "\n"
                elif isinstance(f, GenomeFeature):
                    to_return = to_return + "\tIntergenic Feature:  " + f.protein_accession + "\tLocus Tag: " + f.locus_tag + "\n"
        except:
            to_return = "An error occured printing an operon for " + str(self.genome_fragment_name) +  "(" + str(self.genome_accession) + ")\n"
        
        return to_return