# oprn_consv_calc (Operon Conservation Calculator)

### Overview & How It Works

The `operon_conserve_detect.py` pipeline is a comprehensive tool designed to assess the structural conservation and sequence identity of a reference operon across a specific taxonomic clade. 

**Biological Rationale:**
Operons are clusters of co-transcribed genes that often encode proteins involved in the same biological pathway. Over evolutionary time, operons can undergo rearrangements, gene insertions, deletions, or complete degradation. This pipeline automatically detects putative operon homologs in target genomes and evaluates their conservation by comparing them against a user-defined reference operon.

**Pipeline Algorithm:**
1. **Search & Identification:** The pipeline queries the genes of the reference operon (via protein accessions, locus tags, or raw translations) against a target database using TBLASTN.
2. **Feature Mapping:** It maps the raw BLAST hits to actual biological features (CDS) annotated in the target genomes' GenBank records.
3. **Operon Assembly:** Hits within the same genome fragment are assembled into putative operons based on user-defined proximity thresholds (maximum intergenic distance in base pairs) and feature limits (maximum number of non-homologous intervening genes allowed).
4. **Scoring:** For each assembled operon, the script calculates:
   * **Average Amino Acid Identity (AAI):** A high-performance, lazy-evaluated global alignment score against the reference proteins.
   * **Structural Similarity Score:** A quantitative metric assessing the conservation of gene order, adjacency, and strand orientation relative to the reference operon.
5. **Validation (Reverse BLAST):** To ensure true homology and exclude paralogs, putative hits can be subjected to a reverse BLAST search against the reference genome.

![Workflow](/extra/operon_detect_pipeline.svg)

---

### Citation
If you use this pipeline in your research, please cite the paper where this methodology was first introduced:

> Kitts, G., ... & Erill, I. (2023). The evolutionary dynamics of the LexA regulon in the Vibrionales. *Frontiers in Microbiology, 14*, 1175143. 
> [https://doi.org/10.3389/fmicb.2023.1175143](https://doi.org/10.3389/fmicb.2023.1175143)

*(Note: Additional citations for plasmid-specific analyses will be added upon publication).*

---

### Installation & Dependencies

It is highly recommended to run this pipeline within an isolated `conda` environment. 

#### 1. Python Environment Setup
```bash
# Create a new conda environment
conda create -n oprn_calc python=3.9
conda activate oprn_calc

# Install required Python packages
conda install -c conda-forge biopython pandas tqdm
```

#### 2. Local BLAST+ (Soft Dependency)
The pipeline can operate entirely via the NCBI API (serverless). However, a local installation of the **NCBI BLAST+ command line applications** (`makeblastdb`, `blastp`, `tblastn`) is required **only if**:
* You set `"blast_type": "local"` to search against a local offline database.
* You enable `"reverse_blast": true` (the script builds a temporary local database of your reference genome to verify hit homology).

**Installing BLAST+ via Conda (Recommended):**
```bash
conda install -c bioconda blast
```
Alternatively, you can install BLAST+ directly from the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Ensure that the executables are added to your system's `PATH`.

---

### Usage

Once the environment is set up and the input JSON is configured, run the script via the command line:

```bash
python operon_conserve_detect.py input/your_config_file.json
```

---

### Input Setup

The pipeline is entirely controlled via a single JSON configuration file. See `TEMPLATE_INPUT.json` in the `/input/` directory for a working example.

| Parameter | Description |
|---|---|
| **`entrez` Parameters** | |
| `request_limit` | The max number of attempts that should be made to complete an Entrez request. |
| `sleep_time` | The amount of time that should be waited between consecutive attempts at an Entrez request (seconds). |
| `email` / `api_key` | Email and API key to use for the NCBI Entrez API. |
| **`reference` Parameters** | *(Unified reference configuration)* |
| `database_mode` | Can be `"remote"` or `"local"`. Defines how the reference operon and assembly are loaded. |
| `data` | **Remote mode:** The NCBI nucleotide accession of the reference genome (e.g., `NC_005340.1`). <br>**Local mode:** The file path to the main local GenBank file (`.gbk`) containing the reference operon. |
| `reference_genome_name` | The name of the reference assembly (required for remote mode). |
| **`input_records`** | *(Queries)* |
| `type` | Defines the format of your input queries. Supported types are:<br>• `"protein_accession"`: Standard NCBI protein IDs (e.g., `NP_414878.1`).<br>• `"locus_tag"`: Gene locus tags (e.g., `b345`). Recommended for local `.gbk` files.<br>• `"translation"`: Raw amino acid sequence strings. <br>*Note:* The pipeline utilizes "Input Normalization" to automatically map any input to the reference feature and extract the standard protein accessions. |
| `values` | An array containing the query values in the exact biological order of the reference operon (5' to 3'). |
| **`blast` Parameters** | |
| `tax_include` / `tax_exclude` | NCBI Taxonomy IDs to include or exclude from the BLAST search (e.g., `["562"]` for E. coli). |
| `database` | The database to conduct the BLAST search against (e.g., `"ref_prok_rep_genomes"` or `"nr"`). |
| `blast_type` | Can be `"remote"` or `"local"`. If `"local"`, provide the path prefix in `local_db_path`. |
| `e-val` | The maximum e-value threshold for the BLAST searches. |
| `coverage_min` | The coverage threshold to restrict the BLAST results. |
| `max_hits` | Starting max number of hits to return from a BLAST search. |
| `extensive_search` | If `true`, multiple BLAST searches are conducted iteratively (increasing max hits by `search_mult_factor` up to `max_attempts`) to pull *all* hits within the e-value cutoff. |
| `reverse_blast` | If `true`, all returned hits are tested for true homology by BLAST searching the reference genome with the hit. The local reference database is built automatically in the `./reverse_blast/` directory. |
| **`hit_feature_detection`** | |
| `margin_limit` | The acceptable margin (in bp) between the BLAST alignment positions and the annotated GenBank feature boundaries. |
| `max_attempts` / `mult_factor`| How often and by what factor the `margin_limit` increases if a feature is not immediately found. |
| **`operon_assembly`** | |
| `feature_limit` | Maximum number of non-homologous intervening features allowed between two target genes in an operon. |
| `intergenic_limit` | Maximum distance in bp allowed between two genes in an operon. |
| `use_ref_limit` / `ref_limit_margin` | If `true`, the `intergenic_limit` is dynamically scaled based on the max intergenic distance of the original reference operon multiplied by the margin. |
| **Misc. Parameters** | |
| `thread_limit` | Maximum number of threads for processing species objects. |
| `species_percent_id_limit` | For a species to be outputted, it must have at least one hit with an AAI above this limit. |
| `color_code` | A dictionary defining the HEX colors for each gene in the reference operon used for SVG drawing. `intergenic` refers to inserted non-reference features. |

---

### Output

The pipeline generates a comprehensive, visually rich output folder for each run. By default, results are stored in your defined `output_dir` and include:

**1. The Master CSV (`output.csv`)**
A highly detailed data table containing:
* `Species Name`, `Taxonomic ID`, and `Genome Assembly Accession`.
* The **Structural Similarity Score** (percentage of correctly ordered operon gene pairs).
* The **Average Amino Acid Identity (AAI)** across all genes.
* **Per-Gene Metadata:** For each queried gene, the CSV exports its specific `nucleotide_id`, `strand`, `start`, `stop`, and calculated `protein_accession` to allow for manual biological validation.

**2. Automated iTOL Datasets (`/itol/`)**
The pipeline automatically parses the `output.csv` using NCBI Taxonomy to build taxonomic trees and interactive annotation datasets for the [Interactive Tree Of Life (iTOL)](https://itol.embl.de/). 
*How to use:* Upload the `.nwk` tree to your iTOL account, then simply drag and drop the `.txt` annotation files onto the browser window.
* `taxonomy_tree_full.nwk` / `taxonomy_tree_genus.nwk`: The calculated taxonomic Newick trees.
* `itol_AAI_HEATMAP.txt` & `itol_STRUCT_HEATMAP.txt`: Color gradients for AAI and Structural Similarity.
* `itol_PROTEIN_HEATMAP.txt`: Individual heatmaps for each reference gene.
* `itol_POPUP_INFO.txt`: Interactive HTML popups for each leaf displaying the exact assembly accessions, metadata, and individual percent identities.

**3. SVG Operon Diagrams (`/svg/`)**
A visual genome map for every assembled putative operon across all species. Genes are colored according to the `color_code` dictionary.

---

### Object Definitions (Developer Documentation)

* **GenomeFeature (`features.py`)**
	* Holds information about a biological feature.
	* Member variables: `genome_accession`, `coding_start`, `coding_end`, `strand`, `protein_accession`, `locus_tag`, `aa_sequence`.
	* Main Functions: `get_intergenic_distance()`

* **AnnotatedHit (`features.py`)**
	* A subclass of `GenomeFeature` containing additional BLAST alignment metadata.
	* Member variables: `query_accession`, `align_start`, `align_end`, `percent_identity`, `alignment_seq`.
	* Main Functions: `fetch_feature()` (determines the true GenBank feature corresponding to the BLAST hit).

* **Operon (`operon.py`)**
	* Holds a set of `GenomeFeature`/`AnnotatedHit` objects that have been assembled into the same putative operon.
	* Member variables: `features` (list), `genome_accession`, `strand`.
	* Main Functions: `add_feature()` (adds and dynamically sorts features 5' to 3').

* **GenomeFragment (`genome_fragment.py`)**
	* Represents a single nucleotide accession (e.g., plasmid, chromosome).
	* Member variables: `hits`, `all_features`, `operons`, `genome_accession`, `assembly_accession`, `taxid`.
	* Main Functions: `fetch_features()`, `fetch_record()`, `fetch_hit_features()`, `assemble_operons()`.

* **Species (`species.py`)**
	* Aggregates all `GenomeFragment` objects belonging to the same taxonomic species.
	* Member variables: `assembly_accession`, `species_name`, `genome_fragments`, `sim_score`, `query_percent_ids`.
	* Main Functions: `measure_sim()` (calculates structural similarity), `draw_figure()` (generates SVG diagrams).
