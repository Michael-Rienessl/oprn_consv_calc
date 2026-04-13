'''
@author: ichaudr

'''

from features import AnnotatedHit, GenomeFeature

class Operon:

    def __init__(self, genome_fragment_name, genome_accession, genome_features, strand):
        self.features = []
        self.genome_fragment_name = genome_fragment_name
        self.genome_accession =genome_accession
        self.genome_features = genome_features
        self.strand = strand

        self.local_sim = 0.0
        self.local_aai = 0.0

    def calculate_local_stats(self, reference_pairs):
        '''
        Calculates local structural similarity (SIM)
        and average amino acid identity (AAI) just for this fragment.
        '''

        #1 Calculate local AAI (average of all AnnotedHits in this operon)
        hits = [f for f in self.features if isinstance(f, AnnotatedHit)]
        if hits:
            self.local_aai = sum(float(h.percent_identity) for h in hits) / len(hits)
        else:
            self.local_aai = 0.0

        #2 Calculate local SIM
        # Check how many reference pairs exist in dies specific fragment
        if not reference_pairs:
            self.local_sim = 0.0
            return
        
        found_pairs_count = 0

        # Extract all pairs of query accessions that are directly adjacent in this operon
        # (Intergenic features are ignored because they do not break the chain)
        query_sequence = [f.query_accession for f in self.features if isinstance(f, AnnotatedHit)]

        local_pairs = []
        for i in range(len(query_sequence) - 1):
            local_pairs.append((query_sequence[i], query_sequence[i+1]))
            # Also count reverse pairs (since orientation in the genome can vary)
            local_pairs.append((query_sequence[i+1], query_sequence[i]))

        for ref_p in reference_pairs:
            if ref_p in local_pairs:
                found_pairs_count += 1

        # Local SIM = Found pairs in this fragment / Total reference pairs
        self.local_sim = found_pairs_count / len(reference_pairs)
    
    def add_feature(self, feature):
        '''
        Appends feature to the list of features associated with this operon and sorts it from 5' to 3'

        Parameters
        ----------
        feature: GenomeFeature object
            Feature to be added
        
        Returns
        -------
        None
        '''
        #Checking if the feature is present in the features associated with the genome for this operon.
        if feature in self.genome_features:
            self.features.append(feature)
            self.features = sorted(self.features, key=lambda feature: feature.five_end)
        #else:
         #   raise Exception("The feature you are trying to add is not in the genome assigned for this operon.")
    
    
    
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