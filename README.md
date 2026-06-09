# oprn_consv_calc (Operon Conservation Calculator)

## Overview & How It Works

The `operon_conserve_detect.py` pipeline is a comprehensive tool designed to assess the **structural conservation** and **sequence identity** of a user-defined reference operon across nucleotide databases, genomes, plasmids, contigs, or local GenBank collections.

The pipeline combines:
* Protein-to-nucleotide homology search using **TBLASTN**
* Mapping of raw BLAST hits back to annotated **CDS features** in GenBank records
* Assembly of nearby homologous genes into putative operon fragments
* Structural conservation scoring and Amino-Acid Identity (AAI) scoring
* Optional reverse-BLAST validation
* Automated output tables and SVG operon diagrams
* Automated iTOL-compatible **Taxonomic trees** with rich metadata popups

The pipeline supports both **remote NCBI-based analyses** and **local GenBank / local BLAST database workflows**.

---

## Main Features

* **Flexible Search Modes:** Remote NCBI BLAST or local BLAST+ searches.
* **Multiple Reference Types:** Define operons via remote NCBI nucleotide records, local GenBank files, protein accessions, locus tags, or raw amino-acid translations.
* **Precise Feature Mapping:** Maps BLAST hits back to annotated GenBank CDS features.
* **Structural Assembly:** Assembles putative operons based on strand, gene order, intergenic distance, and allowed intervening CDS features.
* **Permutation Discovery:** Supports strict reference-order scoring or an optional permutation mode for rearranged operons.
* **Automated iTOL Integration:** Connects to NCBI Entrez to download exact taxonomic relationships and generates ready-to-upload iTOL datasets (heatmaps and HTML popups).

---

## Repository Structure

Typical relevant files:

```text
oprn_consv_calc/
├── operon_conserve_detect.py
├── features.py
├── genome_fragment.py
├── species.py
├── operon.py
├── itol_pipeline/
│   └── build_itol_tree_and_datasets.py
├── input/
│   ├── TEMPLATE_INPUT.json
│   └── your_config.json
├── input_data/
│   └── reference_or_target_files.gbk
├── cache/
├── reverse_blast/
└── output/
```

---

## Installation & Dependencies

It is highly recommended to run this pipeline inside an isolated conda environment.

```bash
# Create and activate environment
conda create -n oprn_calc python=3.9
conda activate oprn_calc

# Install Python dependencies
conda install -c conda-forge biopython pandas tqdm

# Install BLAST+
conda install -c bioconda blast
```

The following command-line tools should be available in your `PATH`:
`makeblastdb`, `blastp`, `tblastn`.

---

## Quick Start

Place your JSON configuration file in the `input/` directory (e.g., `input/my_config.json`) and run the pipeline from the project root using only the config file name:

```bash
python operon_conserve_detect.py my_config.json
```

*(Important: The script prepends `input/` internally. Do not type the folder name in the command).*

---

## Input JSON Schema & Parameters

The pipeline uses a unified input schema. Here is an explanation of the most critical configuration blocks.

### 1. Reference & Input Records
Defines where the reference comes from and which genes are queried.

```json
"reference": {
  "database_mode": "remote",
  "data": "NC_005340.1",
  "reference_genome_name": "Rhodopseudomonas_palustris_Test"
},
"input_records": {
  "type": "protein_accession",
  "values": ["NP_958095.1", "NP_958096.1", "NP_958094.1"],
  "names": ["mucC", "mucA", "mucB"]
}
```
* **`database_mode`**: Can be `"remote"` (fetches from NCBI) or `"local"` (uses a local GenBank file path in `"data"`).
* **`values`**: The specific query sequences/accessions.
* **`names`**: Human-readable names used in the output CSVs (e.g., `mucC`, `mucA`, `mucB`).

### 2. Operon Assembly & Permutations
This block controls how the pipeline decides if two genes belong to the same operon.

```json
"operon_assembly": [{
  "feature_limit": 1,
  "intergenic_limit": 150,
  "allow_permutations": true
}]
```
* **`intergenic_limit`**: Maximum allowed distance (in base pairs) between two adjacent genes.
* **`feature_limit`**: Maximum number of intervening, non-query CDS features allowed between two target genes.
* **`allow_permutations` (Crucial Feature)**: 
  * If set to **`false`**: The pipeline strictly enforces the biological order of the reference operon (e.g., `mucC` -> `mucA` -> `mucB`). Any deviation (like `mucA` -> `mucC`) will score lower or be rejected.
  * If set to **`true`**: The pipeline calculates the score for **all possible gene arrangements**. If a plasmid carries a rearranged but intact operon (e.g., `mucA` -> `mucC` -> `mucB`), the pipeline will detect it, score it as highly conserved, and report the discovered order in the output CSV under the `Best_Permutation` column.

---

## Output Files & Data Dictionary

The pipeline generates an organized output directory for each run containing table data and SVG images. The two most important files are the summary and detailed CSVs.

### 1. The Summary Table (`output_summary.csv`)
This table provides a high-level overview. There is exactly **one row per species or assembly**.

| Column | Description |
|---|---|
| `Species` | Organism name extracted from GenBank metadata. |
| `TaxID` | NCBI Taxonomy ID. |
| `Assembly_Accession` | Assembly accession or fallback nucleotide accession. |
| `Best_Operon_Nucleotide_ID` | The specific nucleotide record (e.g., plasmid or chromosome) containing the best-preserved operon fragment. |
| `Reference_Order` | The expected gene order based on the user reference. |
| `Best_Operon_SIM` | Structural Similarity (0-100%) of the *best single fragment* found in this organism (Strict order). |
| `Total_Genomic_SIM` | Global conservation score accounting for all fragments, duplicates, and scattered genes across the entire genome. |
| `Best_Permutation` | *(If allow_permutations=true)* The actual gene order observed in the best fragment (e.g., `mucA, mucC, mucB`). |
| `Best_Permutation_SIM` | *(If allow_permutations=true)* The structural similarity score matching the rearranged order. |
| `Total_Permutation_SIM` | *(If allow_permutations=true)* Global conservation score matching the rearranged order. |
| `Best_Operon_Avg_AAI` | Mean Amino-Acid Identity of the genes within the best fragment. |
| `Global_Total_Avg_AAI` | Mean AAI across all detected hits in the entire genome. |
| `Best_Op_[Gene]_AAI` | Specific AAI for each individual query gene in the best operon. |
| `Total_Count_[Gene]` | Total number of detected paralogs/copies for a specific query gene in the genome. |

### 2. The Detailed Table (`output_detailed.csv`)
A complete, raw export of **every single hit and assembled fragment** detected. Multiple rows per species can exist.

| Column | Description |
|---|---|
| `Species` & `Nucleotide Accession` | Organism and the specific molecule (Plasmid/Chromosome) the fragment is on. |
| `SIM_float` & `Fragment_SIM` | The structural similarity of this specific fragment (e.g., 50.0%). |
| `Fragment_Avg_AAI` | The mean AAI of this specific fragment. |
| `[Gene]_AAI` | The Amino-Acid Identity of the hit. |
| `[Gene]_Accession` | The matched GenBank protein accession or locus tag. |
| `[Gene]_Start` & `[Gene]_Stop` | Exact genomic coordinates of the CDS. |
| `[Gene]_Strand` | Strand orientation (`+` or `-`). |
| `Genetic_Order` | The actual order of genes as they appear on the nucleotide sequence. |
| `Pattern` | A string representing the presence of genes (e.g., `mucA-mucB`). |

*(Note: This table is essential for manual validation of complex rearrangements, partial operons, or scattered gene modules).*

---

## Automated iTOL Tree Generation

To visualize your results, the pipeline includes a standalone script that connects to the NCBI Taxonomy database, builds a robust phylogenetic/taxonomic tree, and generates custom drag-and-drop datasets (Heatmaps and HTML Popups) for [iTOL (Interactive Tree Of Life)](https://itol.embl.de/).

**Usage:**
```bash
python itol_pipeline/build_itol_tree_and_datasets.py path/to/output_summary.csv --mode local --email your.email@example.com
```

This creates an `itol_out/taxonomic/` directory containing:
* `taxonomy_tree_full.nwk` (The Newick Tree file)
* `itol_BEST_SIM_HEATMAP.txt` (Heatmap for Operon intactness)
* `itol_PROTEIN_HEATMAP.txt` (Heatmap for Amino Acid Identities)
* `itol_POPUP_INFO.txt` (Interactive HTML popups containing gene counts).

*(Note: Missing data is automatically handled via iTOL's required `X` syntax to ensure error-free uploads. If iTOL reports "Couldn't find ID..." during upload, this is expected behavior due to NCBI subspecies naming conventions and will not break the visualization).*

---

## Citation

If you use this pipeline in your research, please cite:

> Kitts, G., ... & Erill, I. (2023). *The evolutionary dynamics of the LexA regulon in the Vibrionales.* Frontiers in Microbiology, 14, 1175143.  
> [https://doi.org/10.3389/fmicb.2023.1175143](https://doi.org/10.3389/fmicb.2023.1175143)

Please also cite any external databases and tools used in your analysis, including NCBI BLAST, NCBI Entrez, and iTOL.