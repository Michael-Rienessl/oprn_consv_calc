#!/usr/bin/env python3
"""
Build a taxonomy tree and iTOL annotation datasets from the operon_conserve_detect output table.

This script combines two tasks in one workflow:

1) Build a taxonomy-based Newick tree from the Taxonomic IDs in the output CSV
   using NCBI Entrez Taxonomy.
2) Generate iTOL-compatible dataset files for visualizing:
   - Average Amino Acid Identity (AAI)
   - Operon Structural Similarity
   - Per-reference-protein values
   - Popup information for each leaf

Typical input:
- output.csv from operon_conserve_detect.py

Typical output:
- taxonomy_tree_full.nwk
- taxonomy_tree_genus.nwk
- itol_AAI_HEATMAP.txt
- itol_STRUCT_HEATMAP.txt
- itol_PROTEIN_HEATMAP.txt
- itol_POPUP_INFO.txt

Expected CSV columns:
- Species Name
- Structural Similarity
- Average Percent Amino Acid Identity
- Taxonomic ID
- Genome Assembly Accession
- one or more per-reference-protein columns

Tree leaf labels are written as:
    Species_name__taxid12345

Example:
    python build_itol_tree_and_datasets.py output.csv --email your@email.org --outdir itol_out
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
from io import StringIO

import pandas as pd
from Bio import Entrez, Phylo
from Bio.Phylo.Newick import Tree, Clade


DEFAULT_EDGE_LEN = 1.0
DEFAULT_SLEEP = 0.34
TAXID_RE = re.compile(r"__taxid(\d+)$")


# -------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------

def clean_name(s: str) -> str:
    s = str(s).strip()
    s = s.replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", s)
    return s


def pct_to_float(value) -> float:
    """Convert values like '84.2%' or '84.2' to float. Returns NaN if conversion fails."""
    if pd.isna(value):
        return float("nan")
    s = str(value).strip()
    if s == "":
        return float("nan")
    if s.lower() in {"none", "na", "n/a", "nan", "null"}:
        return float("nan")
    if s.endswith("%"):
        s = s[:-1]
        
    try:
        return float(s)
    except ValueError:
        return float("nan")


def extract_leaf_labels(newick_text: str) -> List[str]:
    tree = Phylo.read(StringIO(newick_text), "newick")
    return [clade.name.strip() for clade in tree.get_terminals() if clade.name]


def build_taxid_to_leafid(labels: List[str]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for lbl in labels:
        match = TAXID_RE.search(lbl)
        if match:
            taxid = int(match.group(1))
            mapping.setdefault(taxid, lbl)
    return mapping


# -------------------------------------------------------------------
# Tree-building code
# -------------------------------------------------------------------

def read_taxids(csv_path: Path) -> List[int]:
    taxids = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "Taxonomic ID" not in reader.fieldnames:
            raise KeyError(f"CSV has no column 'Taxonomic ID'. Columns: {reader.fieldnames}")
        for row in reader:
            v = str(row["Taxonomic ID"]).strip()
            if v:
                taxids.append(int(v))
    return sorted(set(taxids))


def fetch_taxonomy_records(taxids: List[int], sleep_s: float) -> list:
    records_all = []
    batch_size = 200

    for i in range(0, len(taxids), batch_size):
        batch = taxids[i:i + batch_size]
        handle = Entrez.efetch(
            db="taxonomy",
            id=",".join(map(str, batch)),
            retmode="xml"
        )
        recs = Entrez.read(handle)
        handle.close()
        records_all.extend(recs)
        time.sleep(sleep_s)

    return records_all


def build_tax_graph(tax_records: list, leaf_taxids: List[int]):
    """
    Build parent->children using LineageEx.
    Also capture taxid->name and taxid->rank.
    """
    parent = {}
    children = defaultdict(set)
    name = {}
    rank = {}

    for rec in tax_records:
        tid = int(rec["TaxId"])
        name[tid] = clean_name(rec.get("ScientificName", str(tid)))
        rank[tid] = rec.get("Rank", "no_rank")

        lineage = rec.get("LineageEx", [])
        prev = None
        for node in lineage:
            ntid = int(node["TaxId"])
            name[ntid] = clean_name(node.get("ScientificName", str(ntid)))
            rank.setdefault(ntid, node.get("Rank", "no_rank"))

            if prev is not None:
                parent[ntid] = prev
                children[prev].add(ntid)
            prev = ntid

        if lineage:
            last = int(lineage[-1]["TaxId"])
            parent[tid] = last
            children[last].add(tid)

    used_nodes = set(name.keys())
    roots = [n for n in used_nodes if n not in parent]

    if len(roots) == 1:
        root_id = roots[0]
    else:
        root_id = -1
        name[root_id] = "ROOT"
        rank[root_id] = "no_rank"
        for r in roots:
            parent[r] = root_id
            children[root_id].add(r)

    for tid in leaf_taxids:
        name.setdefault(tid, str(tid))
        rank.setdefault(tid, "no_rank")

    return root_id, children, name, rank


def to_newick_clade(node_id, children, name, leaf_taxids_set, edge_len=DEFAULT_EDGE_LEN):
    """
    Recursively build clades and assign a positive branch_length to every edge.
    """
    ch = sorted(children.get(node_id, []))
    if not ch:
        label = name.get(node_id, str(node_id))
        if node_id in leaf_taxids_set:
            label = f"{label}__taxid{node_id}"
        return Clade(name=label)

    clade_children = []
    for c in ch:
        child_clade = to_newick_clade(c, children, name, leaf_taxids_set, edge_len=edge_len)
        child_clade.branch_length = float(edge_len)
        clade_children.append(child_clade)

    label = name.get(node_id, str(node_id))
    return Clade(clades=clade_children, name=label)


def collapse_to_genus(root_id, children, rank):
    """
    Build a reduced child-map where we stop expanding below genus.
    """
    new_children = defaultdict(set)

    def rec(n):
        for c in children.get(n, []):
            if rank.get(c) == "genus":
                new_children[n].add(c)
            else:
                new_children[n].add(c)
                rec(c)

    rec(root_id)

    pruned = defaultdict(set)

    def rec2(n):
        for c in new_children.get(n, []):
            pruned[n].add(c)
            if rank.get(c) != "genus":
                rec2(c)

    rec2(root_id)
    return pruned


def build_taxonomy_trees(csv_path: Path, outdir: Path, sleep_s: float) -> Tuple[Path, Path]:
    leaf_taxids = read_taxids(csv_path)
    print(f"Loaded {len(leaf_taxids)} unique TaxIDs.")

    tax_recs = fetch_taxonomy_records(leaf_taxids, sleep_s=sleep_s)
    print(f"Fetched {len(tax_recs)} taxonomy records.")

    root_id, children, name, rank = build_tax_graph(tax_recs, leaf_taxids)
    leaf_set = set(leaf_taxids)

    out_full_nwk = outdir / "taxonomy_tree_full.nwk"
    out_genus_nwk = outdir / "taxonomy_tree_genus.nwk"

    root_clade = to_newick_clade(root_id, children, name, leaf_set, edge_len=DEFAULT_EDGE_LEN)
    full_tree = Tree(root=root_clade)
    Phylo.write(full_tree, str(out_full_nwk), "newick")

    genus_children = collapse_to_genus(root_id, children, rank)
    root_clade_g = to_newick_clade(root_id, genus_children, name, leaf_set, edge_len=DEFAULT_EDGE_LEN)
    genus_tree = Tree(root=root_clade_g)
    Phylo.write(genus_tree, str(out_genus_nwk), "newick")

    print(f"Wrote: {out_full_nwk}")
    print(f"Wrote: {out_genus_nwk}")

    return out_full_nwk, out_genus_nwk


# -------------------------------------------------------------------
# iTOL dataset generation
# -------------------------------------------------------------------

def detect_metric_columns(df: pd.DataFrame) -> Tuple[str, str, List[str]]:
    """
    Detect the core metrics and per-protein columns.
    """
    aai_candidates = [
        "Average Percent Amino Acid Identity",
        "AAI",
        "AAI_mean",
        "Average Amino Acid Identity",
    ]
    struct_candidates = [
        "Structural Similarity",
        "Operon Structural Similarity",
        "Operon structural similarity",
        "Structure Similarity",
        "Structural_similarity",
    ]

    aai_col = next((c for c in aai_candidates if c in df.columns), None)
    struct_col = next((c for c in struct_candidates if c in df.columns), None)

    if aai_col is None:
        raise KeyError("Could not find an AAI column in the CSV.")
    if struct_col is None:
        raise KeyError("Could not find a structural similarity column in the CSV.")

    metadata_cols = {
        "Species Name",
        "Taxonomic ID",
        "Genome Assembly Accession",
        "Sequence Accessions",
        "taxid",
        "leaf_id",
        aai_col,
        struct_col,
    }

    protein_cols = [c for c in df.columns if c not in metadata_cols]

    valid_protein_cols = []
    for col in protein_cols:
        converted = df[col].apply(pct_to_float)
        if converted.notna().any():
            valid_protein_cols.append(col)

    return aai_col, struct_col, valid_protein_cols


def write_heatmap_header(
    handle,
    dataset_label: str,
    color: str,
    field_labels: List[str],
    min_value: float = 0,
    max_value: float = 100,
    color_min: str = "#ff0000",
    color_max: str = "#00ff00",
):
    handle.write("DATASET_HEATMAP\n")
    handle.write("SEPARATOR TAB\n")
    handle.write(f"DATASET_LABEL\t{dataset_label}\n")
    handle.write(f"COLOR\t{color}\n")
    handle.write("SHOW_INTERNAL\t0\n")
    handle.write("AUTO_LEGEND\t1\n")
    handle.write(f"FIELD_LABELS\t" + "\t".join(field_labels) + "\n")
    handle.write(f"COLOR_MIN\t{color_min}\n")
    handle.write(f"COLOR_MAX\t{color_max}\n")
    handle.write(f"USER_MIN_VALUE\t{min_value}\n")
    handle.write(f"USER_MAX_VALUE\t{max_value}\n")
    handle.write("DATA\n")


def write_aai_heatmap(df_ok: pd.DataFrame, outpath: Path, aai_col: str):
    with outpath.open("w", encoding="utf-8") as f:
        write_heatmap_header(
            f,
            dataset_label="AAI (%)",
            color="#1f4db3",
            field_labels=["AAI"],
            min_value=0,
            max_value=100,
        )
        for _, row in df_ok.iterrows():
            aai = pct_to_float(row[aai_col])
            if pd.notna(aai):
                f.write(f"{row['leaf_id']}\t{aai:.4f}\n")


def write_struct_heatmap(df_ok: pd.DataFrame, outpath: Path, struct_col: str):
    with outpath.open("w", encoding="utf-8") as f:
        write_heatmap_header(
            f,
            dataset_label="Operon Structural Similarity (%)",
            color="#b012e6",
            field_labels=["Structure"],
            min_value=0,
            max_value=100,
        )
        for _, row in df_ok.iterrows():
            struct = pct_to_float(row[struct_col])
            if pd.notna(struct):
                f.write(f"{row['leaf_id']}\t{struct:.4f}\n")


def write_protein_heatmap(df_ok: pd.DataFrame, outpath: Path, protein_cols: List[str]):
    if not protein_cols:
        return

    with outpath.open("w", encoding="utf-8") as f:
        write_heatmap_header(
            f,
            dataset_label="Reference Protein Similarities (%)",
            color="#4c4c4c",
            field_labels=protein_cols,
            min_value=0,
            max_value=100,
        )
        for _, row in df_ok.iterrows():
            values = []
            for col in protein_cols:
                val = pct_to_float(row[col])
                if pd.isna(val):
                    values.append("")
                else:
                    values.append(f"{val:.4f}")
            f.write(f"{row['leaf_id']}\t" + "\t".join(values) + "\n")


def write_popup_info(
    df_ok: pd.DataFrame,
    outpath: Path,
    aai_col: str,
    struct_col: str,
    protein_cols: List[str],
):
    with outpath.open("w", encoding="utf-8") as f:
        f.write("POPUP_INFO\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATA\n")

        for _, row in df_ok.iterrows():
            species = str(row.get("Species Name", "NA"))
            taxid = row.get("Taxonomic ID", "NA")
            assembly = str(row.get("Genome Assembly Accession", "NA"))

            aai = pct_to_float(row[aai_col])
            struct = pct_to_float(row[struct_col])

            lines = [
                f"<b>{html.escape(species)}</b>",
                f"TaxID: {html.escape(str(taxid))}",
                f"Assembly: {html.escape(assembly)}",
                f"AAI: {aai:.2f}%" if pd.notna(aai) else "AAI: NA",
                f"Operon Structural Similarity: {struct:.2f}%" if pd.notna(struct) else "Operon Structural Similarity: NA",
            ]

            if protein_cols:
                lines.append("<br><b>Reference protein values</b>")
                for col in protein_cols:
                    val = pct_to_float(row[col])
                    if pd.notna(val):
                        lines.append(f"{html.escape(col)}: {val:.2f}%")
                    else:
                        lines.append(f"{html.escape(col)}: NA")

            popup_html = "<br>".join(lines)
            title = html.escape(species)

            f.write(f"{row['leaf_id']}\t{title}\t{popup_html}\n")


def load_and_match(tree_path: Path, csv_path: Path) -> Tuple[pd.DataFrame, List[str]]:
    newick = tree_path.read_text(encoding="utf-8")
    labels = extract_leaf_labels(newick)
    taxid_to_leaf = build_taxid_to_leafid(labels)

    if not taxid_to_leaf:
        raise ValueError(
            "Could not extract any TaxIDs from tree leaf labels. "
            "Expected labels like 'Species_name__taxid12345'."
        )

    df = pd.read_csv(csv_path)

    if "Taxonomic ID" not in df.columns:
        raise KeyError("CSV must contain the column 'Taxonomic ID'.")
    if "Species Name" not in df.columns:
        raise KeyError("CSV must contain the column 'Species Name'.")

    df["taxid"] = df["Taxonomic ID"].astype(int)
    df["leaf_id"] = df["taxid"].map(taxid_to_leaf)

    return df, labels


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build a taxonomy tree and iTOL datasets from operon_conserve_detect output.csv"
    )
    parser.add_argument("csv", help="Path to output.csv")
    parser.add_argument("--outdir", default="itol_out", help="Output directory")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", default=None, help="Optional NCBI API key")
    parser.add_argument(
        "--sleep",
        type=float,
        default=DEFAULT_SLEEP,
        help=f"Delay between Entrez batch requests (default: {DEFAULT_SLEEP})",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path.resolve()}")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # 1) Build trees
    full_tree_path, genus_tree_path = build_taxonomy_trees(
        csv_path=csv_path,
        outdir=outdir,
        sleep_s=args.sleep,
    )

    # 2) Build iTOL datasets from full tree
    df, labels = load_and_match(full_tree_path, csv_path)
    aai_col, struct_col, protein_cols = detect_metric_columns(df)

    missing = df[df["leaf_id"].isna()][["Species Name", "taxid"]]
    if len(missing) > 0:
        print("WARNING: Some TaxIDs from the CSV were not found in the tree (showing up to 20 rows):")
        print(missing.head(20).to_string(index=False))

    df_ok = df.dropna(subset=["leaf_id"]).copy()

    aai_path = outdir / "itol_AAI_HEATMAP.txt"
    struct_path = outdir / "itol_STRUCT_HEATMAP.txt"
    protein_path = outdir / "itol_PROTEIN_HEATMAP.txt"
    popup_path = outdir / "itol_POPUP_INFO.txt"

    write_aai_heatmap(df_ok, aai_path, aai_col)
    write_struct_heatmap(df_ok, struct_path, struct_col)
    write_protein_heatmap(df_ok, protein_path, protein_cols)
    write_popup_info(df_ok, popup_path, aai_col, struct_col, protein_cols)

    print(f"Wrote: {aai_path}")
    print(f"Wrote: {struct_path}")
    if protein_cols:
        print(f"Wrote: {protein_path}")
        print("Protein columns included in heatmap:")
        for col in protein_cols:
            print(f"  - {col}")
    else:
        print("No per-protein columns detected; protein heatmap was not populated.")
    print(f"Wrote: {popup_path}")
    print(f"Matched {len(df_ok)}/{len(df)} CSV rows to tree leaves.")
    print(f"Tree leaves detected: {len(labels)}")
    print(f"Full tree: {full_tree_path}")
    print(f"Genus tree: {genus_tree_path}")


if __name__ == "__main__":
    main()