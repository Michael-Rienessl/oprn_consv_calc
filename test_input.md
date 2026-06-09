# `input.json` Parameter Reference

This document explains the parameters used by the `oprn_consv_calc` input JSON file. The pipeline reads the JSON in `operon_conserve_detect.py` via `load_input_file()`, then uses the values to define the reference operon, BLAST search mode, hit annotation, operon assembly, output generation, and NCBI Entrez behavior.

The input file is JSON. Most top-level entries are dictionaries or one-element lists of dictionaries, for example:

```json
"blast": [{
  "blast_type": "local",
  "local_db_path": "./plsdb_plasmids_database/my_nucl_db"
}]
```

---

## 1. Minimal high-level structure

A typical input file contains these top-level blocks:

```json
{
  "reference": {},
  "input_records": {},
  "entrez": [{}],
  "blast": [{}],
  "hit_feature_detection": [{}],
  "operon_assembly": [{}],
  "thread_limit": 4,
  "species_percent_id_limit": 0.45,
  "output_dir": "output/my_run/",
  "cache_dir": "cache/my_run/",
  "outputs": {},
  "color_code": [{}]
}
```

---

## 2. Supported analysis scenarios

The pipeline separates the **reference source** from the **target BLAST database**.

| Scenario | `reference.database_mode` | `blast.blast_type` | Meaning |
|---|---|---|---|
| Remote reference vs remote BLAST | `"remote"` | `"remote"` | Reference genome/proteins are obtained from NCBI; target search uses NCBI BLAST. |
| Remote reference vs local BLAST | `"remote"` | `"local"` | Reference is obtained from NCBI; target search uses a local nucleotide BLAST database. |
| Local reference vs remote BLAST | `"local"` | `"remote"` | Reference genes are extracted from a local GenBank file; target search uses NCBI BLAST. |
| Local reference vs local BLAST | `"local"` | `"local"` | Reference genes are extracted from a local GenBank file; target search uses a local nucleotide BLAST database. |

---

# `reference`

The `reference` block defines the genome or GenBank file from which the reference operon is defined.

```json
"reference": {
  "database_mode": "local",
  "data": "input_data/SI_3H_plasmid.gbk",
  "reference_genome_name": ""
}
```

## `reference.database_mode`

- Type: `string`
- Allowed values: `"local"`, `"remote"`
- Required: yes

Determines where the reference genome record comes from.

### `"local"`

The reference record is read from a local GenBank file.

```json
"database_mode": "local",
"data": "input_data/SI_3H_plasmid.gbk"
```

Use this when your reference operon is in a local `.gb`, `.gbk`, `.gbff`, or `.genbank` file.

### `"remote"`

The reference record is fetched from NCBI using a nucleotide accession.

```json
"database_mode": "remote",
"data": "NC_005340.1",
"reference_genome_name": "Rhodopseudomonas_palustris_Test"
```

Use this when the reference genome or plasmid is available from NCBI.

---

## `reference.data`

- Type: `string`
- Required: yes

The content depends on `reference.database_mode`.

| Mode | Expected value |
|---|---|
| `"local"` | Path to a local GenBank file |
| `"remote"` | NCBI nucleotide accession, including version if possible |

Examples:

```json
"data": "input_data/SI_3H_plasmid.gbk"
```

```json
"data": "NC_005340.1"
```

---

## `reference.reference_genome_name`

- Type: `string`
- Required in `"remote"` reference mode
- Optional / can be empty in `"local"` reference mode

A user-defined name for the reference genome. In remote mode, the current parser requires this field.

Example:

```json
"reference_genome_name": "Rhodopseudomonas_palustris_Test"
```

---

## `reference.reference_genome_assembly`

- Type: `string`
- Required: no
- Default: `""`

Optional NCBI assembly accession for the remote reference genome.

If this is omitted in remote mode, the pipeline attempts to resolve the assembly from the reference nucleotide accession.

Example:

```json
"reference_genome_assembly": "GCF_000000000.0"
```

---

# `input_records`

The `input_records` block defines the genes/proteins that form the reference operon.

```json
"input_records": {
  "type": "locus_tag",
  "values": [
    "AJMIMBIC_04697",
    "AJMIMBIC_04698",
    "AJMIMBIC_04699"
  ],
  "names": [
    "mucC",
    "mucA",
    "mucB"
  ]
}
```

## `input_records.type`

- Type: `string`
- Allowed values used by the current pipeline: `"locus_tag"`, `"protein_accession"`, `"protein_id"`, `"translation"`; helper logic also supports `"gene"` in some paths.
- Required: yes

Defines how the entries in `input_records.values` should be interpreted.

Recommended usage:

| Reference mode | Recommended `input_records.type` | Example |
|---|---|---|
| Local GenBank reference | `"locus_tag"` | `"AJMIMBIC_04697"` |
| Remote NCBI reference | `"protein_id"` / `"protein_accession"` | `"NP_958095.1"` |
| Raw amino-acid input | `"translation"` | `"MKKLL..."` |

Notes:

- In local reference mode, `"locus_tag"` is usually the safest option because local GenBank annotations often contain stable locus tags.
- For remote protein queries, use NCBI protein accessions such as `NP_...`, `WP_...`, or `YP_...`.
- The code contains support for both `"protein_id"` and `"protein_accession"` patterns in different helper sections. For consistency, use `"protein_id"` when your JSON follows the examples using NCBI protein IDs.

---

## `input_records.values`

- Type: `list[string]`
- Required: yes
- Minimum length: `2`

The ordered list of reference operon genes/proteins.

Example:

```json
"values": [
  "AJMIMBIC_04697",
  "AJMIMBIC_04698",
  "AJMIMBIC_04699"
]
```

The order matters because the pipeline calculates structural similarity from adjacent gene pairs.

For example, with:

```text
mucC - mucA - mucB
```

the expected adjacent pairs are:

```text
mucC-mucA
mucA-mucB
```

---

## `input_records.names`

- Type: `list[string]`
- Required: no
- Default: each value is used as its own name

Human-readable names for output tables and plots.

Example:

```json
"names": [
  "mucC",
  "mucA",
  "mucB"
]
```

If `names` is provided, the output columns use these names, for example:

```text
Best_Op_mucC_AAI
Total_Count_mucC
```

instead of accession/locus-tag-based names.

---

# `entrez`

The `entrez` block controls NCBI Entrez requests.

```json
"entrez": [{
  "request_limit": 3,
  "sleep_time": 1,
  "email": "your.email@example.com",
  "api_key": "YOUR_NCBI_API_KEY"
}]
```

## `entrez[0].request_limit`

- Type: `integer`
- Allowed values: positive integers
- Recommended: `3`
- Required: yes

Maximum number of retry attempts for NCBI requests.

Used when downloading GenBank records, protein FASTA sequences, assembly metadata, and other Entrez resources.

---

## `entrez[0].sleep_time`

- Type: `float` or `integer`
- Allowed values: `>= 0`
- Recommended: `0.34` to `1`
- Required: yes

Pause in seconds between NCBI requests or retry attempts.

For remote-heavy analyses, higher values reduce the risk of NCBI throttling or HTTP 429 errors.

---

## `entrez[0].email`

- Type: `string`
- Required: yes for NCBI usage

Email address passed to NCBI Entrez.

Use a real email address for NCBI compliance.

Example:

```json
"email": "your.email@example.com"
```

---

## `entrez[0].api_key`

- Type: `string`
- Required: no, but strongly recommended for remote analyses

NCBI API key used by Entrez.

Do not commit real API keys to public repositories.

Example:

```json
"api_key": "YOUR_NCBI_API_KEY"
```

---

# `blast`

The `blast` block controls the target search.

```json
"blast": [{
  "blast_type": "local",
  "local_genomes_dir": "./PLSDB/PLSDB/gbk",
  "local_db_path": "./plsdb_plasmids_database/my_nucl_db",
  "database": "ref_prok_rep_genomes",
  "e_val": 1e-10,
  "coverage_min": 0.7,
  "max_hits": 1000000,
  "max_attempts": 3,
  "search_mult_factor": 2,
  "annotate": true,
  "extensive_search": true,
  "reverse_blast": true
}]
```

## `blast[0].blast_type`

- Type: `string`
- Allowed values: `"remote"`, `"local"`
- Required: yes
- Default in code: `"remote"`

Defines whether target search uses NCBI BLAST or a local BLAST database.

---

## `blast[0].database`

- Type: `string`
- Required for remote BLAST
- Default in code: `"ref_prok_rep_genomes"`

NCBI BLAST database name for remote searches.

Common examples:

```json
"database": "ref_prok_rep_genomes"
```

```json
"database": "nt"
```

For local BLAST, this value is not the active search database; `local_db_path` is used instead.

---

## `blast[0].local_db_path`

- Type: `string`
- Required when `blast_type` is `"local"`
- Ignored for remote BLAST

Path prefix of the local nucleotide BLAST database.

Important: this is the BLAST database prefix, not a single file.

Example:

```json
"local_db_path": "./plsdb_plasmids_database/my_nucl_db"
```

This expects files such as:

```text
my_nucl_db.nhr
my_nucl_db.nin
my_nucl_db.nsq
```

or equivalent BLAST database files.

---

## `blast[0].local_genomes_dir`

- Type: `string`
- Required for local BLAST if GenBank annotation is needed
- Required if the pipeline should auto-build the local BLAST database
- Ignored for remote BLAST

Folder containing the local GenBank records corresponding to the local BLAST database.

Example:

```json
"local_genomes_dir": "./PLSDB/PLSDB/gbk"
```

The GenBank record IDs should match the local BLAST FASTA headers.

In strict local mode, the intended relationship is:

```text
BLAST hit ID == FASTA header == GenBank record.id
```

---

## `blast[0].tax_include`

- Type: `list[string]` or `list[int]`
- Used by: remote BLAST
- Ignored by: local BLAST unless the local database was pre-filtered manually
- Default: `[]`

Restrict remote BLAST to one or more NCBI taxonomy IDs.

Example:

```json
"tax_include": ["1236"]
```

This restricts NCBI BLAST to Gammaproteobacteria.

For Enterobacteriaceae:

```json
"tax_include": ["543"]
```

---

## `blast[0].tax_exclude`

- Type: `list[string]` or `list[int]`
- Used by: remote BLAST
- Default: `[]`

Exclude one or more NCBI taxonomy IDs from remote BLAST.

Example:

```json
"tax_exclude": ["562"]
```

---

## `blast[0].e_val`

- Type: `float`
- Recommended: `1e-10`
- Required: yes

Maximum BLAST E-value for accepted hits.

Example:

```json
"e_val": 1e-10
```

Lower values are stricter.

---

## `blast[0].coverage_min`

- Type: `float`
- Allowed values: usually `0.0` to `1.0`
- Required: yes
- Recommended:
  - `0.7` to `0.85` for strict full-length homolog detection
  - `0.3` to `0.5` for exploratory or divergent searches

Minimum fraction of the query protein covered by the BLAST HSP.

Example:

```json
"coverage_min": 0.7
```

---

## `blast[0].max_hits`

- Type: `integer`
- Required: yes

Maximum number of target sequences returned by BLAST.

For remote BLAST, this becomes the NCBI BLAST `hitlist_size`.

For local BLAST, this is used as `max_target_seqs`.

Examples:

```json
"max_hits": 2000
```

```json
"max_hits": 1000000
```

For large local plasmid databases, a high value may be needed to avoid truncating hits.

---

## `blast[0].max_attempts`

- Type: `integer`
- Used mainly by: remote BLAST
- Default in code: `3`

Maximum number of repeated BLAST attempts.

Remote BLAST may increase `max_hits` if the result appears incomplete.

---

## `blast[0].search_mult_factor`

- Type: `integer` or `float`
- Used mainly by: remote BLAST
- Default in code: `2`

If a remote BLAST search appears incomplete, `max_hits` is multiplied by this factor for the next attempt.

Example:

```json
"search_mult_factor": 2
```

---

## `blast[0].annotate`

- Type: `boolean`
- Required: should be `true`
- Default in code: `true`

If `true`, BLAST hits are stored as `AnnotatedHit` objects.

Important: downstream grouping and operon assembly currently require annotated hits. Keep this as:

```json
"annotate": true
```

---

## `blast[0].extensive_search`

- Type: `boolean`
- Used mainly by: remote BLAST
- Default in code: `true`

If `true`, remote BLAST can retry/increase hit limits to improve completeness.

If `false`, the search is only attempted once.

---

## `blast[0].reverse_blast`

- Type: `boolean`
- Default in code: `false`

If `true`, each candidate hit is checked by reverse BLAST against a protein database built from the reference genome.

This improves specificity but can strongly increase runtime.

Recommended for final analyses:

```json
"reverse_blast": true
```

Recommended for quick tests:

```json
"reverse_blast": false
```

---

## `blast[0].allow_permutations`

- Type: `boolean`
- Current recommended location: `operon_assembly[0].allow_permutations`

Some older or experimental JSONs place this field under `blast`, but the current parser uses `allow_permutations` from the `operon_assembly` block.

Use:

```json
"operon_assembly": [{
  "allow_permutations": true
}]
```

---

# `hit_feature_detection`

This block controls how BLAST HSPs are mapped to annotated CDS features.

```json
"hit_feature_detection": [{
  "margin_limit": 5,
  "max_attempts": 3,
  "mult_factor": 3
}]
```

## `hit_feature_detection[0].margin_limit`

- Type: `integer`
- Default examples: `5`, `15`
- Required: yes

Historical parameter for feature matching tolerance.

In the current SeqRecord-based feature matching, the primary logic uses overlap/coverage between the BLAST hit coordinates and CDS coordinates. The field is still parsed and passed through for compatibility.

---

## `hit_feature_detection[0].max_attempts`

- Type: `integer`
- Required: yes

Historical retry parameter for feature assignment.

Still parsed and passed through for compatibility.

---

## `hit_feature_detection[0].mult_factor`

- Type: `integer` or `float`
- Required: yes

Historical multiplier used to expand feature-search tolerance in older feature assignment logic.

Still parsed and passed through for compatibility.

---

# `operon_assembly`

This block controls how hits are grouped into operon fragments and how structural conservation is scored.

```json
"operon_assembly": [{
  "feature_limit": 1,
  "intergenic_limit": 150,
  "use_ref_limit": true,
  "ref_limit_margin": 0.5,
  "allow_permutations": true
}]
```

## `operon_assembly[0].feature_limit`

- Type: `integer`
- Required: yes
- Recommended:
  - `0` for no intervening features
  - `1` for strict but tolerant operon assembly
  - `3` for more permissive assembly

Maximum number of non-query CDS features that may be inserted between adjacent query hits while still considering them part of the same operon fragment.

Example:

```json
"feature_limit": 1
```

---

## `operon_assembly[0].intergenic_limit`

- Type: `integer`
- Unit: base pairs
- Required: yes

Maximum allowed distance between neighboring features for them to be connected into one operon fragment.

Example:

```json
"intergenic_limit": 150
```

---

## `operon_assembly[0].use_ref_limit`

- Type: `boolean`
- Required: yes

If `true`, the pipeline calculates the maximum intergenic distance in the reference operon and uses it to set the working `intergenic_limit`.

The derived value is:

```text
intergenic_limit = (1 + ref_limit_margin) × max_reference_intergenic_distance
```

---

## `operon_assembly[0].ref_limit_margin`

- Type: `float`
- Required if `use_ref_limit` is `true`
- Recommended: `0.5`

Margin added to the reference operon intergenic distance when `use_ref_limit` is active.

Example:

```json
"ref_limit_margin": 0.5
```

means 50% extra tolerance relative to the largest reference intergenic gap.

---

## `operon_assembly[0].allow_permutations`

- Type: `boolean`
- Default in code: `false`
- Recommended:
  - `false` for strict reference-order conservation
  - `true` for exploratory discovery of alternative gene orders

If `false`, structural similarity is calculated only against the input/reference gene order.

If `true`, the pipeline evaluates all permutations of the reference gene order and reports the best permutation-based score.

This is useful when the same genes are present but the order may be rearranged.

Outputs affected when enabled:

```text
Best_Permutation
Best_Permutation_SIM
Total_Permutation_SIM
```

---

# `thread_limit`

```json
"thread_limit": 4
```

- Type: `integer`
- Default in code: `1`

Intended maximum number of processing threads.

Important implementation note: in the current lazy-loading species-processing block, this value is parsed but not fully enforced as a global thread pool limit. Fragment threads are created per assembly/species group. Therefore, this parameter is currently more of a compatibility/configuration field than a strict global CPU limiter.

---

# `species_percent_id_limit`

```json
"species_percent_id_limit": 0.45
```

- Type: `float`
- Allowed values: usually `0.0` to `1.0`
- Recommended: `0.45`

Filtering threshold applied after feature extraction and identity recalculation.

A species/assembly/plasmid group is retained if at least one query gene has a maximum identity greater than or equal to this threshold.

Example:

```json
"species_percent_id_limit": 0.45
```

means that a group passes if at least one reference gene has a hit with >=45% amino-acid identity.

---

# `output_dir`

```json
"output_dir": "output/plsdb_plasmid_analyses_final_new/"
```

- Type: `string`
- Required: recommended
- Default in code: `"output/"`

Run-specific output folder.

The pipeline writes outputs such as:

```text
output_summary.csv
output_detailed.csv
pipeline_run.log
svg/
itol/
```

depending on configuration.

The run log is written into this folder as:

```text
pipeline_run.log
```

---

# `cache_dir`

```json
"cache_dir": "./remote_cache/gbk"
```

- Type: `string`
- Required: recommended
- Default in code: `"./cache/"`

Folder used for cached GenBank records and other intermediate files.

Behavior depends on the scenario:

| Scenario | Recommended `cache_dir` |
|---|---|
| Remote BLAST | A dedicated writable cache folder, e.g. `"./remote_cache/gbk"` |
| Local BLAST with local GBKs | Can point to the local GenBank folder or a separate cache |
| Local reference only | Any writable cache folder |

For remote downloads, cached GenBank records are expected to use accession-based filenames such as:

```text
CP000000.1.gb
```

---

# `outputs`

This optional block controls expensive downstream output generation.

```json
"outputs": {
  "generate_svg_graphs": false,
  "generate_taxonomic_trees": false,
  "generate_phylogenetic_trees": false
}
```

## `outputs.generate_svg_graphs`

- Type: `boolean`
- Default in code: `true`
- Recommended for large runs: `false`

Controls whether SVG operon diagrams are generated.

If `true`, the pipeline creates SVG diagrams under:

```text
output_dir/svg/
```

For very large runs this can be slow and generate many files.

---

## `outputs.generate_taxonomic_trees`

- Type: `boolean`
- Default in code: `true`
- Recommended for large or test runs: `false`

Controls whether taxonomic iTOL tree/dataset generation is run.

If `false`, the iTOL/taxonomic tree pipeline is skipped.

---

## `outputs.generate_phylogenetic_trees`

- Type: `boolean`
- Default in code: `false`

Controls whether Sourmash-based phylogenetic tree generation is run.

The Sourmash code is currently kept in the pipeline but disabled by default.

Use:

```json
"generate_phylogenetic_trees": true
```

only if the Sourmash workflow and input files are ready.

---

# `color_code`

The `color_code` block defines colors for SVG operon diagrams.

```json
"color_code": [{
  "intergenic": "#000000",
  "AJMIMBIC_04697": "#ffa200",
  "AJMIMBIC_04698": "#68ff03",
  "AJMIMBIC_04699": "#05eeff"
}]
```

## `color_code[0].intergenic`

- Type: hex color string
- Required if SVG graphs are enabled

Color used for intervening non-query CDS features in SVG diagrams.

Example:

```json
"intergenic": "#000000"
```

---

## Gene-specific color entries

- Type: hex color string
- Required if SVG graphs are enabled
- Key should match the query accession/locus tag used as `feat.query_accession`

Each query gene should have a color.

Example:

```json
"AJMIMBIC_04697": "#ffa200"
```

If SVG generation is disabled, missing or incomplete color definitions usually do not matter.

---

# Example 1: Local reference vs remote BLAST

```json
{
  "reference": {
    "database_mode": "local",
    "data": "input_data/SI_3H_plasmid.gbk",
    "reference_genome_name": ""
  },

  "input_records": {
    "type": "locus_tag",
    "values": [
      "AJMIMBIC_04697",
      "AJMIMBIC_04698",
      "AJMIMBIC_04699"
    ],
    "names": [
      "mucC",
      "mucA",
      "mucB"
    ]
  },

  "entrez": [{
    "request_limit": 3,
    "sleep_time": 1,
    "email": "your.email@example.com",
    "api_key": "YOUR_NCBI_API_KEY"
  }],

  "blast": [{
    "blast_type": "remote",
    "tax_include": ["1236"],
    "tax_exclude": [],
    "database": "ref_prok_rep_genomes",
    "e_val": 1e-10,
    "coverage_min": 0.85,
    "max_hits": 50000,
    "max_attempts": 3,
    "search_mult_factor": 2,
    "annotate": true,
    "extensive_search": true,
    "reverse_blast": true
  }],

  "hit_feature_detection": [{
    "margin_limit": 5,
    "max_attempts": 3,
    "mult_factor": 3
  }],

  "operon_assembly": [{
    "feature_limit": 1,
    "intergenic_limit": 150,
    "use_ref_limit": true,
    "ref_limit_margin": 0.5,
    "allow_permutations": true
  }],

  "thread_limit": 4,
  "species_percent_id_limit": 0.45,
  "output_dir": "output/local_reference_vs_remote_blast/",
  "cache_dir": "./remote_cache/gbk",

  "outputs": {
    "generate_svg_graphs": false,
    "generate_taxonomic_trees": false,
    "generate_phylogenetic_trees": false
  },

  "color_code": [{
    "intergenic": "#000000",
    "AJMIMBIC_04697": "#ffa200",
    "AJMIMBIC_04698": "#68ff03",
    "AJMIMBIC_04699": "#05eeff"
  }]
}
```

---

# Example 2: Local reference vs local BLAST database

```json
{
  "reference": {
    "database_mode": "local",
    "data": "input_data/SI_3H_plasmid.gbk",
    "reference_genome_name": ""
  },

  "input_records": {
    "type": "locus_tag",
    "values": [
      "AJMIMBIC_04697",
      "AJMIMBIC_04698",
      "AJMIMBIC_04699"
    ],
    "names": [
      "mucC",
      "mucA",
      "mucB"
    ]
  },

  "entrez": [{
    "request_limit": 3,
    "sleep_time": 1,
    "email": "your.email@example.com",
    "api_key": "YOUR_NCBI_API_KEY"
  }],

  "blast": [{
    "blast_type": "local",
    "local_genomes_dir": "./PLSDB/PLSDB/gbk",
    "local_db_path": "./plsdb_plasmids_database/my_nucl_db",
    "database": "ref_prok_rep_genomes",
    "e_val": 1e-10,
    "coverage_min": 0.7,
    "max_hits": 1000000,
    "max_attempts": 3,
    "search_mult_factor": 2,
    "annotate": true,
    "extensive_search": true,
    "reverse_blast": true
  }],

  "hit_feature_detection": [{
    "margin_limit": 5,
    "max_attempts": 3,
    "mult_factor": 3
  }],

  "operon_assembly": [{
    "feature_limit": 1,
    "intergenic_limit": 150,
    "use_ref_limit": true,
    "ref_limit_margin": 0.5,
    "allow_permutations": true
  }],

  "thread_limit": 4,
  "species_percent_id_limit": 0.45,
  "output_dir": "output/local_reference_vs_local_blast/",
  "cache_dir": "./PLSDB/PLSDB/gbk",

  "outputs": {
    "generate_svg_graphs": false,
    "generate_taxonomic_trees": false,
    "generate_phylogenetic_trees": false
  },

  "color_code": [{
    "intergenic": "#000000",
    "AJMIMBIC_04697": "#ffa200",
    "AJMIMBIC_04698": "#68ff03",
    "AJMIMBIC_04699": "#05eeff"
  }]
}
```

---

# Example 3: Remote reference vs remote BLAST

```json
{
  "reference": {
    "database_mode": "remote",
    "data": "NC_005340.1",
    "reference_genome_name": "Rhodopseudomonas_palustris_Test"
  },

  "input_records": {
    "type": "protein_id",
    "values": [
      "NP_958095.1",
      "NP_958096.1",
      "NP_958094.1"
    ],
    "names": [
      "gene1",
      "gene2",
      "gene3"
    ]
  },

  "entrez": [{
    "request_limit": 3,
    "sleep_time": 1,
    "email": "your.email@example.com",
    "api_key": "YOUR_NCBI_API_KEY"
  }],

  "blast": [{
    "tax_include": ["543"],
    "tax_exclude": [],
    "database": "ref_prok_rep_genomes",
    "blast_type": "remote",
    "local_db_path": "",
    "e_val": 1e-10,
    "coverage_min": 0.3,
    "max_hits": 2000,
    "max_attempts": 3,
    "search_mult_factor": 2,
    "annotate": true,
    "extensive_search": true,
    "reverse_blast": true
  }],

  "hit_feature_detection": [{
    "margin_limit": 5,
    "max_attempts": 3,
    "mult_factor": 3
  }],

  "operon_assembly": [{
    "feature_limit": 1,
    "intergenic_limit": 150,
    "use_ref_limit": true,
    "ref_limit_margin": 0.5,
    "allow_permutations": true
  }],

  "thread_limit": 4,
  "species_percent_id_limit": 0.45,
  "output_dir": "output/remote_reference_vs_remote_blast/",
  "cache_dir": "cache_remote_vs_remote/",

  "outputs": {
    "generate_svg_graphs": false,
    "generate_taxonomic_trees": false,
    "generate_phylogenetic_trees": false
  },

  "color_code": [{
    "intergenic": "#000000",
    "NP_958095.1": "#ff0000",
    "NP_958094.1": "#ff0002",
    "NP_958096.1": "#ff0001"
  }]
}
```

---

# Example 4: Remote reference vs local BLAST database

```json
{
  "reference": {
    "database_mode": "remote",
    "data": "NC_005340.1",
    "reference_genome_name": "Rhodopseudomonas_palustris_Test"
  },

  "input_records": {
    "type": "protein_id",
    "values": [
      "NP_958095.1",
      "NP_958096.1",
      "NP_958094.1"
    ],
    "names": [
      "gene1",
      "gene2",
      "gene3"
    ]
  },

  "entrez": [{
    "request_limit": 3,
    "sleep_time": 1,
    "email": "your.email@example.com",
    "api_key": "YOUR_NCBI_API_KEY"
  }],

  "blast": [{
    "blast_type": "local",
    "local_genomes_dir": "./PLSDB/PLSDB/gbk",
    "local_db_path": "./plsdb_plasmids_database/my_nucl_db",
    "database": "ref_prok_rep_genomes",
    "e_val": 1e-10,
    "coverage_min": 0.5,
    "max_hits": 1000000,
    "max_attempts": 3,
    "search_mult_factor": 2,
    "annotate": true,
    "extensive_search": true,
    "reverse_blast": true
  }],

  "hit_feature_detection": [{
    "margin_limit": 5,
    "max_attempts": 3,
    "mult_factor": 3
  }],

  "operon_assembly": [{
    "feature_limit": 1,
    "intergenic_limit": 150,
    "use_ref_limit": true,
    "ref_limit_margin": 0.5,
    "allow_permutations": true
  }],

  "thread_limit": 4,
  "species_percent_id_limit": 0.45,
  "output_dir": "output/remote_reference_vs_local_blast/",
  "cache_dir": "./PLSDB/PLSDB/gbk",

  "outputs": {
    "generate_svg_graphs": false,
    "generate_taxonomic_trees": false,
    "generate_phylogenetic_trees": false
  },

  "color_code": [{
    "intergenic": "#000000",
    "NP_958095.1": "#ff0000",
    "NP_958094.1": "#ff0002",
    "NP_958096.1": "#ff0001"
  }]
}
```

---

# Practical recommendations

## For quick tests

Use:

```json
"max_hits": 100,
"reverse_blast": false,
"outputs": {
  "generate_svg_graphs": false,
  "generate_taxonomic_trees": false,
  "generate_phylogenetic_trees": false
}
```

## For final local PLSDB-scale runs

Use:

```json
"blast_type": "local",
"max_hits": 1000000,
"reverse_blast": true,
"outputs": {
  "generate_svg_graphs": false,
  "generate_taxonomic_trees": false,
  "generate_phylogenetic_trees": false
}
```

Enable SVGs only after the main CSV output is correct.

## For remote BLAST runs

Use smaller `max_hits` first, then increase if needed.

Remote BLAST can be slow due to NCBI queueing, network latency, Entrez metadata resolution, and reverse BLAST.

---

# Output files

The main output files are:

## `output_summary.csv`

One row per retained assembly/species/plasmid group, depending on the run mode and metadata grouping.

Important columns include:

```text
Species
TaxID
Assembly_Accession
Best_Operon_Nucleotide_ID
Reference_Order
Best_Operon_SIM
Total_Genomic_SIM
Best_Operon_Avg_AAI
Global_Total_Avg_AAI
Best_Op_<gene>_AAI
Total_Count_<gene>
```

If `allow_permutations` is enabled, additional columns are written:

```text
Best_Permutation
Best_Permutation_SIM
Total_Permutation_SIM
```

## `output_detailed.csv`

Fragment-level table containing every assembled operon fragment.

Important columns include:

```text
Species
TaxID
Assembly_Accession
Nucleotide Accession
Fragment_SIM
Fragment_Avg_AAI
<gene>_AAI
<gene>_Accession
<gene>_Start
<gene>_Stop
<gene>_Strand
Genetic_Order
```

## `pipeline_run.log`

Run-specific log file written into `output_dir`.

Example:

```text
output/my_run/pipeline_run.log
```

## Optional output folders

Depending on the `outputs` block:

```text
svg/
itol/
```
