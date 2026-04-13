# oprn_consv_calc (Operon Conservation Calculator)

## Overview & How It Works
The `operon_conserve_detect.py` pipeline is a comprehensive tool designed to assess the **structural conservation** and **sequence identity** of a reference operon across a specific taxonomic clade.

---

## Biological Rationale
Operons are clusters of co-transcribed genes that often encode proteins involved in the same biological pathway. Over evolutionary time, operons can undergo rearrangements, gene insertions, deletions, or complete duplication.

This pipeline automatically detects putative operon homologs and evaluates their conservation by comparing them against a user-defined reference, providing a distinction between:

- the **best-preserved operon** (single best fragment per genome), and
- the **overall genomic presence** of target genes (global score that accounts for paralogs and noise).

---

## Pipeline Algorithm

### 1) Search & Identification
Queries reference genes via **TBLASTN** against target databases.

### 2) Feature Mapping
Maps raw BLAST hits to biological features (**CDS**) in GenBank records.

### 3) Operon Assembly
Assembles hits into putative operons based on:
- proximity thresholds, and
- intervening gene limits.

### 4) Scoring (Dual-Layer Analysis)

- **Best Operon SIM / AAI**  
  Identifies the *single most conserved operon fragment* in the genome.

- **Total Genomic SIM**  
  Calculates a global architecture score normalized by gene frequency to account for paralogs and genomic noise.

### 5) Memory Management (Lazy Writing)
To support large-scale analyses, species data is exported to a detailed log and then purged from RAM immediately after processing (`clean()` function).

---

## Citation
If you use this pipeline in your research, please cite:

> Kitts, G., ... & Erill, I. (2023). *The evolutionary dynamics of the LexA regulon in the Vibrionales.* Frontiers in Microbiology, 14, 1175143.  
> https://doi.org/10.3389/fmicb.2023.1175143

---

## Installation & Dependencies
It is highly recommended to run this pipeline within an isolated conda environment.

```bash
# Create and activate environment
conda create -n oprn_calc python=3.9
conda activate oprn_calc

# Install dependencies
conda install -c conda-forge biopython pandas tqdm
conda install -c bioconda blast
```

---

## Usage
Once the environment is set up and the input JSON is configured, run the script via:

```bash
python operon_conserve_detect.py input/your_config_file.json
```

---

## Output Structure
The pipeline generates an organized output directory for each run:

### 1) The Summary Table (`output_summary.csv`)
The primary file for comparative analysis, listing one row per species:

- **Best_Operon_SIM** & **Best_Operon_Avg_AAI**: quality metrics for the single highest-scoring operon fragment.
- **Total_Genomic_SIM** & **Global_Total_Avg_AAI**: global conservation metrics reflecting the entire genome.
- **Best_Op_[Gene]_AAI**: individual identity scores for every gene specifically within the "Best Operon".
- **Total_Count_[Gene]**: genomic copy number (counts) for each reference gene across the whole genome.

### 2) The Detailed Table (`output_detailed.csv`)
A complete, “on-the-fly” export of every single hit and fragment detected. Essential for manual validation of fragmented operons or complex gene duplications.

### 3) Automated iTOL Datasets (`/itol/`)
- `taxonomy_tree_full.nwk`: Newick tree of all detected species.
- `itol_BEST_SIM_HEATMAP.txt`: visualization of the highest quality operon per species.
- `itol_TOTAL_SIM_HEATMAP.txt`: visualization of global genomic architecture conservation.
- `itol_PROTEIN_HEATMAP.txt`: multi-column heatmap for individual gene identities (Best Operon).
- `itol_POPUP_INFO.txt`: interactive popups with assembly IDs, metadata, and genomic gene counts.

### 4) SVG Operon Diagrams (`/svg/`)
Visual genome maps of every assembled operon fragment, color-coded by gene identity.

---

## Object Definitions (Developer Documentation)

- **Operon** (`operon.py`)  
  Calculates `local_sim` and `local_aai` for specific fragments using a neighbor-pair comparison algorithm.

- **Species** (`species.py`)  
  Implements the global scoring logic (`measure_sim`) and critical memory cleanup (`clean()`).

- **GenomeFragment** (`genome_fragment.py`)  
  Handles the assembly of BLAST hits into logical operon units.

- **AnnotatedHit** (`features.py`)  
  Subclass of `GenomeFeature` containing BLAST alignment metadata and identity scores.