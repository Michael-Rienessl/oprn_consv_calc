#!/usr/bin/env python3

import argparse
import csv
import html
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import Entrez, Phylo
from Bio.Phylo.Newick import Tree, Clade


# ==========================================
# GLOBAL SETTINGS
# ==========================================

DEFAULT_SLEEP = 0.34


# ==========================================
# BASIC HELPERS
# ==========================================

def clean_name(s: str) -> str:
    """
    iTOL/Newick-safe label cleaning.
    Used mostly for phylogenetic/local tree datasets.
    Taxonomic tree IDs use taxid-only labels instead.
    """
    s = str(s).strip().replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9_\-\.]+", "_", s)


def pretty_species_name(s: str) -> str:
    """
    Display-friendly species name for iTOL labels.
    """
    if pd.isna(s):
        return ""

    s = str(s).strip()

    if s == "" or s.lower() in {"nan", "none", "na", "n/a", "null"}:
        return ""

    # If aggregated species strings exist, use first one for the visible label.
    s = s.split(";")[0].strip()

    return s.replace("_", " ")


def detect_separator(path: Path) -> str:
    """
    Auto-detect CSV separator from first line.
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return ";" if first.count(";") > first.count(",") else ","


def pct_to_float(value) -> float:
    """
    Converts percentage strings like '100.0%' or numeric values to float.
    Returns NaN for missing / invalid values.
    """
    if pd.isna(value):
        return float("nan")

    s = str(value).strip()

    if s == "" or s.lower() in {"none", "na", "n/a", "nan", "null", "x"}:
        return float("nan")

    if s.endswith("%"):
        s = s[:-1]

    try:
        return float(s)
    except ValueError:
        return float("nan")


def safe_int_taxid(value) -> Optional[int]:
    """
    Converts TaxID values robustly.
    Returns None for invalid/missing TaxIDs.
    """
    if pd.isna(value):
        return None

    s = str(value).strip()

    if s == "" or s.lower() in {"none", "na", "n/a", "nan", "null"}:
        return None

    try:
        taxid = int(float(s))
    except Exception:
        return None

    if taxid <= 0:
        return None

    return taxid


def taxid_leaf_label(taxid: int) -> str:
    """
    Stable leaf ID used both in taxonomy tree and taxonomy datasets.
    """
    return f"taxid{int(taxid)}"


# ==========================================
# TAXONOMY FUNCTIONS
# ==========================================

def read_taxids(csv_path: Path) -> List[int]:
    """
    Reads unique valid TaxIDs from input CSV.
    """
    sep = detect_separator(csv_path)
    df = pd.read_csv(csv_path, sep=sep)

    if "TaxID" not in df.columns:
        return []

    taxids = []

    for v in df["TaxID"]:
        tid = safe_int_taxid(v)
        if tid is not None:
            taxids.append(tid)

    return sorted(set(taxids))


def fetch_taxonomy_records(taxids: List[int], sleep_s: float) -> list:
    """
    Fetches NCBI taxonomy records in chunks.
    """
    records_all = []

    for i in range(0, len(taxids), 200):
        chunk = taxids[i:i + 200]
        print(f"  Fetching TaxIDs {i + 1}-{i + len(chunk)} / {len(taxids)}")

        handle = Entrez.efetch(
            db="taxonomy",
            id=",".join(map(str, chunk)),
            retmode="xml"
        )
        records_all.extend(Entrez.read(handle))
        handle.close()
        time.sleep(sleep_s)

    return records_all


def build_tax_graph(tax_records, requested_leaf_taxids):
    """
    Builds a taxonomy graph from NCBI taxonomy records.

    Important:
    Every requested TaxID gets a forced terminal child named taxid12345.
    This guarantees that every dataset ID exists as a terminal leaf in iTOL.
    """
    parent = {}
    children = defaultdict(set)
    name = {}

    for rec in tax_records:
        tid = int(rec["TaxId"])
        name[tid] = clean_name(rec.get("ScientificName", str(tid)))

        prev = None

        for node in rec.get("LineageEx", []):
            ntid = int(node["TaxId"])
            name[ntid] = clean_name(node.get("ScientificName", str(ntid)))

            if prev is not None:
                parent[ntid] = prev
                children[prev].add(ntid)

            prev = ntid

        if rec.get("LineageEx", []):
            last = int(rec["LineageEx"][-1]["TaxId"])
            parent[tid] = last
            children[last].add(tid)

    roots = [n for n in name.keys() if n not in parent]

    if len(roots) == 1:
        root_id = roots[0]
    else:
        root_id = -1
        name[root_id] = "ROOT"
        for r in roots:
            children[root_id].add(r)

    for tid in requested_leaf_taxids:
        name.setdefault(tid, str(tid))

    # Critical iTOL fix:
    # Force every requested TaxID to appear as a terminal leaf with a stable taxid-only label.
    # This avoids mismatches when the TaxID itself becomes an internal node.
    for tid in requested_leaf_taxids:
        forced_leaf = f"leaf_taxid{tid}"
        children[tid].add(forced_leaf)
        name[forced_leaf] = taxid_leaf_label(tid)

    return root_id, children, name


def to_newick_clade(node_id, children, name):
    """
    Converts taxonomy graph to a Newick clade.

    Leaves are named exactly as stored in name[node_id].
    For forced leaves this is taxid12345.
    """
    ch = sorted(children.get(node_id, []), key=lambda x: str(x))

    if not ch:
        return Clade(name=name.get(node_id, str(node_id)))

    clade_children = [to_newick_clade(c, children, name) for c in ch]

    for c in clade_children:
        c.branch_length = 1.0

    return Clade(clades=clade_children, name=name.get(node_id, str(node_id)))


def build_taxonomy_tree(csv_path: Path, outdir: Path, sleep_s: float):
    """
    Builds taxonomy tree from all unique TaxIDs in the CSV.
    """
    leaf_taxids = read_taxids(csv_path)

    if not leaf_taxids:
        print("[WARN] No valid TaxIDs found. Skipping taxonomy tree.")
        return None

    print(f"Fetching Taxonomy for {len(leaf_taxids)} unique TaxIDs...")
    tax_recs = fetch_taxonomy_records(leaf_taxids, sleep_s)

    root_id, children, name = build_tax_graph(tax_recs, leaf_taxids)

    out_nwk = outdir / "taxonomy_tree_full.nwk"
    tree = Tree(root=to_newick_clade(root_id, children, name))
    Phylo.write(tree, str(out_nwk), "newick")

    print(f"[OK] Taxonomy tree written: {out_nwk}")
    return out_nwk


# ==========================================
# DATASET COLUMN DETECTION
# ==========================================

def detect_metric_columns(df: pd.DataFrame):
    """
    Detects main metric columns.
    """
    best_sim_col = "Best_Operon_SIM" if "Best_Operon_SIM" in df.columns else None
    total_sim_col = "Total_Genomic_SIM" if "Total_Genomic_SIM" in df.columns else None

    best_op_protein_cols = [
        c for c in df.columns
        if c.startswith("Best_Op_") and c.endswith("_AAI")
    ]

    return best_sim_col, total_sim_col, best_op_protein_cols


def merge_unique_nonempty(values, max_items=8):
    """
    Compact metadata merger for duplicated leaf IDs.
    """
    cleaned = []

    for v in values:
        if pd.isna(v):
            continue

        s = str(v).strip()

        if s == "" or s.lower() in {"nan", "none", "na", "n/a", "null"}:
            continue

        if s not in cleaned:
            cleaned.append(s)

    if not cleaned:
        return ""

    if len(cleaned) > max_items:
        return "; ".join(cleaned[:max_items]) + f"; ... +{len(cleaned) - max_items} more"

    return "; ".join(cleaned)


def aggregate_by_leaf_id(df: pd.DataFrame, best_col, tot_col, prot_cols):
    """
    Ensures one row per leaf_id.

    Taxonomic tree:
    Many plasmids can map to the same TaxID, so we collapse them.

    Aggregation strategy:
    - SIM / AAI metrics: max
    - gene counts: sum
    - metadata strings: unique compact list
    - N_Records: number of source rows collapsed into this leaf
    """
    if "leaf_id" not in df.columns:
        raise KeyError("df must contain leaf_id before aggregation.")

    count_cols = [c for c in df.columns if c.startswith("Total_Count_")]
    metadata_cols = [
        c for c in [
            "Species",
            "TaxID",
            "Mobility",
            "Deepomics_Mobility",
            "Incompatibility",
            "Incompatibility_Type",
            "ESBL",
            "LexA-Box",
            "LexA-Motif",
            "Other_LexA_Target_Count",
            "Other_LexA_Target_Products",
            "Extra_PolV_Count",
            "Extra_PolV_Detail",
        ]
        if c in df.columns
    ]

    rows = []

    for leaf_id, grp in df.groupby("leaf_id", dropna=False):
        row = {"leaf_id": leaf_id, "N_Records": len(grp)}

        for col in [best_col, tot_col] + list(prot_cols):
            if col and col in grp.columns:
                vals = grp[col].apply(pct_to_float)
                max_val = vals.max(skipna=True)
                row[col] = max_val if pd.notna(max_val) else float("nan")

        for col in count_cols:
            vals = pd.to_numeric(grp[col], errors="coerce")
            row[col] = int(vals.sum(skipna=True)) if vals.notna().any() else ""

        for col in metadata_cols:
            row[col] = merge_unique_nonempty(grp[col].tolist())

        rows.append(row)

    return pd.DataFrame(rows)


# ==========================================
# iTOL DATASET WRITERS
# ==========================================

def write_heatmap_header(handle, label: str, color: str, fields: list, vmin=0, vmax=100):
    handle.write("DATASET_HEATMAP\n")
    handle.write("SEPARATOR TAB\n")
    handle.write(f"DATASET_LABEL\t{label}\n")
    handle.write(f"COLOR\t{color}\n")
    handle.write("SHOW_INTERNAL\t0\n")
    handle.write("AUTO_LEGEND\t1\n")
    handle.write("FIELD_LABELS\t" + "\t".join(fields) + "\n")
    handle.write("COLOR_MIN\t#ff0000\n")
    handle.write("COLOR_MAX\t#00ff00\n")
    handle.write(f"USER_MIN_VALUE\t{vmin}\n")
    handle.write(f"USER_MAX_VALUE\t{vmax}\n")
    handle.write("DATA\n")


def write_taxon_label_dataset(df_tax: pd.DataFrame, outdir: Path, include_taxid_in_label: bool = False):
    """
    Writes an iTOL LABELS file.

    Internal tree leaf:
        taxid562

    Displayed label:
        Escherichia coli

    This keeps iTOL matching stable while displaying species names.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "itol_TAXON_LABELS.txt"

    with out.open("w", encoding="utf-8") as f:
        f.write("LABELS\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATA\n")

        for _, r in df_tax.iterrows():
            leaf_id = str(r["leaf_id"])
            species = pretty_species_name(r.get("Species", ""))
            taxid = str(r.get("TaxID", "")).split(";")[0].strip()

            if species:
                label = species
            else:
                label = leaf_id

            if include_taxid_in_label and taxid:
                label = f"{label} ({taxid})"

            f.write(f"{leaf_id}\t{label}\n")

    print(f"[OK] Taxon label dataset saved: {out}")


def write_simple_binary_dataset(df_ok: pd.DataFrame, outdir: Path):
    """
    Optional binary presence track:
    complete / partial / single based on Best_Operon_SIM.
    """
    if "Best_Operon_SIM" not in df_ok.columns:
        return

    out = outdir / "itol_OPERON_CLASS_BINARY.txt"

    with out.open("w", encoding="utf-8") as f:
        f.write("DATASET_BINARY\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tOperon class\n")
        f.write("COLOR\t#333333\n")
        f.write("FIELD_SHAPES\t1\t1\t1\n")
        f.write("FIELD_LABELS\tSingle_or_0SIM\tPartial_50SIM\tComplete_100SIM\n")
        f.write("FIELD_COLORS\t#d95f02\t#7570b3\t#1b9e77\n")
        f.write("DATA\n")

        for _, r in df_ok.iterrows():
            sim = pct_to_float(r.get("Best_Operon_SIM", float("nan")))

            is_0 = 1 if pd.notna(sim) and sim == 0.0 else 0
            is_50 = 1 if pd.notna(sim) and sim == 50.0 else 0
            is_100 = 1 if pd.notna(sim) and sim == 100.0 else 0

            f.write(f"{r['leaf_id']}\t{is_0}\t{is_50}\t{is_100}\n")


def write_datasets(df_ok: pd.DataFrame, outdir: Path, best_col: str, tot_col: str, prot_cols: list):
    """
    Writes iTOL dataset files for the given mapped dataframe.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Best SIM heatmap
    if best_col and best_col in df_ok.columns:
        with (outdir / "itol_BEST_SIM_HEATMAP.txt").open("w", encoding="utf-8") as f:
            write_heatmap_header(f, "Best Operon SIM (%)", "#1f4db3", ["Best SIM"])

            for _, r in df_ok.iterrows():
                v = pct_to_float(r[best_col])
                if pd.notna(v):
                    f.write(f"{r['leaf_id']}\t{v:.4f}\n")

    # Total SIM heatmap
    if tot_col and tot_col in df_ok.columns:
        with (outdir / "itol_TOTAL_SIM_HEATMAP.txt").open("w", encoding="utf-8") as f:
            write_heatmap_header(f, "Total Genomic SIM (%)", "#b012e6", ["Total SIM"])

            for _, r in df_ok.iterrows():
                v = pct_to_float(r[tot_col])
                if pd.notna(v):
                    f.write(f"{r['leaf_id']}\t{v:.4f}\n")

    # Protein AAI heatmap
    if prot_cols:
        with (outdir / "itol_PROTEIN_HEATMAP.txt").open("w", encoding="utf-8") as f:
            write_heatmap_header(f, "Protein Similarities (%)", "#4c4c4c", prot_cols)

            for _, r in df_ok.iterrows():
                vals = []

                for c in prot_cols:
                    val = pct_to_float(r.get(c, float("nan")))

                    if pd.notna(val):
                        vals.append(f"{val:.4f}")
                    else:
                        vals.append("X")

                f.write(f"{r['leaf_id']}\t" + "\t".join(vals) + "\n")

    # Optional binary operon class
    write_simple_binary_dataset(df_ok, outdir)

    # Popup Info
    count_cols = [c for c in df_ok.columns if c.startswith("Total_Count_")]

    extra_cols = [
        "N_Records",
        "Mobility",
        "Deepomics_Mobility",
        "Incompatibility",
        "Incompatibility_Type",
        "ESBL",
        "LexA-Box",
        "LexA-Motif",
        "Other_LexA_Target_Count",
        "Other_LexA_Target_Products",
        "Extra_PolV_Count",
        "Extra_PolV_Detail",
    ]

    found_extras = [c for c in extra_cols if c in df_ok.columns]

    with (outdir / "itol_POPUP_INFO.txt").open("w", encoding="utf-8") as f:
        f.write("POPUP_INFO\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATA\n")

        for _, r in df_ok.iterrows():
            sp = str(r.get("Species", "NA"))
            tax = str(r.get("TaxID", "NA"))

            lines = [
                f"<b>{html.escape(pretty_species_name(sp) or sp)}</b>",
                f"TaxID: {html.escape(tax)}",
                f"Leaf ID: {html.escape(str(r['leaf_id']))}",
            ]

            if "N_Records" in r.index:
                lines.append(f"Collapsed records: {html.escape(str(r.get('N_Records')))}")

            if best_col and best_col in r.index:
                b_sim = pct_to_float(r[best_col])
                lines.append(
                    f"Best Operon SIM: {b_sim:.2f}%"
                    if pd.notna(b_sim)
                    else "Best Operon SIM: NA"
                )

            if tot_col and tot_col in r.index:
                t_sim = pct_to_float(r[tot_col])
                lines.append(
                    f"Total Genomic SIM: {t_sim:.2f}%"
                    if pd.notna(t_sim)
                    else "Total Genomic SIM: NA"
                )

            if prot_cols:
                lines.append("<br><b>Protein Identities:</b>")

                for c in prot_cols:
                    v = pct_to_float(r.get(c, float("nan")))
                    n = c.replace("Best_Op_", "").replace("_AAI", "")

                    if pd.notna(v):
                        lines.append(f"{html.escape(n)}: {v:.2f}%")
                    else:
                        lines.append(f"{html.escape(n)}: NA")

            if count_cols:
                lines.append("<br><b>Gene Counts:</b>")

                for c in count_cols:
                    label = c.replace("Total_Count_", "")
                    lines.append(f"{html.escape(label)}: {html.escape(str(r.get(c, 'NA')))}")

            if found_extras:
                lines.append("<br><b>Metadata:</b>")

                for c in found_extras:
                    v = r.get(c, "")
                    s = str(v)

                    if pd.notna(v) and s not in ["", "nan", "Unknown", "UNKNOWN", "NOT_FOUND"]:
                        lines.append(f"{html.escape(c)}: {html.escape(s)}")

            title = pretty_species_name(sp) or str(r["leaf_id"])
            f.write(f"{r['leaf_id']}\t{html.escape(title)}\t{'<br>'.join(lines)}\n")


# ==========================================
# MAPPINGS
# ==========================================

def map_taxonomic(df: pd.DataFrame, best_col, tot_col, prot_cols):
    """
    Maps rows to taxonomic tree leaf IDs.

    Critical fix:
    leaf_id is taxid-only, e.g. taxid562.
    This exactly matches the forced terminal leaves in taxonomy_tree_full.nwk.
    """
    df_tax = df.copy()

    def make_leaf_id(row):
        tid = safe_int_taxid(row.get("TaxID"))
        if tid is None:
            return None
        return taxid_leaf_label(tid)

    df_tax["leaf_id"] = df_tax.apply(make_leaf_id, axis=1)
    df_tax = df_tax.dropna(subset=["leaf_id"])

    # One taxonomic leaf per TaxID, so collapse duplicate plasmid rows.
    df_tax = aggregate_by_leaf_id(df_tax, best_col, tot_col, prot_cols)

    return df_tax


def map_phylogenetic(df: pd.DataFrame, mode: str):
    """
    Maps rows to phylogenetic tree IDs.

    local mode:
      Species__Best_Operon_Nucleotide_ID_without_version

    remote mode:
      Species

    Use phylogenetic datasets only with the matching phylogenetic tree,
    not with taxonomy_tree_full.nwk.
    """
    df_phylo = df.copy()
    leaf_ids = []

    for _, row in df_phylo.iterrows():
        sp = clean_name(row.get("Species", "Unknown"))

        if mode == "local":
            acc = str(row.get("Best_Operon_Nucleotide_ID", "")).strip()

            if not acc or acc.lower() == "nan":
                acc = str(row.get("Assembly_Accession", "")).strip()

            acc = clean_name(acc.split(".")[0])
            leaf_ids.append(f"{sp}__{acc}")

        else:
            leaf_ids.append(sp)

    df_phylo["leaf_id"] = leaf_ids
    df_phylo = df_phylo.dropna(subset=["leaf_id"])

    if mode == "remote":
        best_col, tot_col, prot_cols = detect_metric_columns(df_phylo)
        df_phylo = aggregate_by_leaf_id(df_phylo, best_col, tot_col, prot_cols)
    else:
        df_phylo = df_phylo.drop_duplicates(subset=["leaf_id"], keep="first")

    return df_phylo


# ==========================================
# VALIDATION
# ==========================================

def read_tree_leaf_ids(tree_path: Path) -> set:
    tree = Phylo.read(str(tree_path), "newick")
    return {term.name for term in tree.get_terminals() if term.name}


def read_dataset_ids(dataset_path: Path) -> set:
    """
    Reads first-column IDs from iTOL dataset DATA section.
    Works for HEATMAP, BINARY, POPUP_INFO and LABELS.
    """
    ids = set()
    in_data = False

    with dataset_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.strip() == "DATA":
                in_data = True
                continue

            if in_data:
                ids.add(line.split("\t")[0])

    return ids


def validate_dataset_ids_against_tree(tree_path: Path, dataset_dir: Path):
    """
    Checks whether all IDs in dataset files exist as terminal leaves in tree.
    """
    if not tree_path or not tree_path.exists():
        print("[WARN] No tree found for validation.")
        return

    tree_ids = read_tree_leaf_ids(tree_path)

    dataset_files = sorted(dataset_dir.glob("itol_*.txt"))

    if not dataset_files:
        print("[WARN] No iTOL dataset files found for validation.")
        return

    print(f"\n--- Validating datasets against tree: {tree_path.name} ---")
    print(f"Tree terminal leaves: {len(tree_ids)}")

    total_missing = 0

    for ds in dataset_files:
        ds_ids = read_dataset_ids(ds)
        missing = sorted(ds_ids - tree_ids)

        print(f"{ds.name}: dataset IDs={len(ds_ids)}, missing={len(missing)}")

        if missing:
            total_missing += len(missing)
            print("  First missing IDs:")
            for x in missing[:20]:
                print(f"    {x}")

    if total_missing == 0:
        print("[OK] All dataset IDs are present in the tree.")
    else:
        print(f"[WARN] Total missing ID occurrences across datasets: {total_missing}")


# ==========================================
# MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(
        description="Build iTOL taxonomy tree and matching dataset files."
    )

    parser.add_argument(
        "csv",
        help="Path to output_summary.csv or output_summary_analyzed.csv"
    )

    parser.add_argument(
        "--outdir",
        default="itol_out",
        help="Output directory"
    )

    parser.add_argument(
        "--email",
        default="",
        help="Email for NCBI Entrez"
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=["local", "remote"],
        help=(
            "local: phylogenetic dataset IDs are Species__Accession. "
            "remote: phylogenetic dataset IDs are Species. "
            "This does not affect the taxonomic tree."
        )
    )

    parser.add_argument(
        "--skip-taxonomy",
        action="store_true",
        help="Skip NCBI taxonomy tree generation."
    )

    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate generated taxonomy datasets against taxonomy_tree_full.nwk."
    )

    parser.add_argument(
        "--label-taxid",
        action="store_true",
        help="Show TaxID in visible iTOL leaf labels, e.g. Escherichia coli (562)."
    )

    args = parser.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)

    tax_dir = outdir / "taxonomic"
    phy_dir = outdir / "phylogenetic"

    tax_dir.mkdir(parents=True, exist_ok=True)
    phy_dir.mkdir(parents=True, exist_ok=True)

    if args.email:
        Entrez.email = args.email
    else:
        print("[WARN] No Entrez email provided. NCBI may reject requests.")
        print("       Use: --email your.name@example.com")

    sep = detect_separator(csv_path)
    df = pd.read_csv(csv_path, sep=sep)

    if "TaxID" not in df.columns:
        raise KeyError("CSV missing required column: TaxID")

    if "Species" not in df.columns:
        raise KeyError("CSV missing required column: Species")

    best_col, tot_col, prot_cols = detect_metric_columns(df)

    print("\n=== Input summary ===")
    print(f"CSV: {csv_path}")
    print(f"Rows: {len(df)}")
    print(f"Separator: {repr(sep)}")
    print(f"Best SIM column: {best_col}")
    print(f"Total SIM column: {tot_col}")
    print(f"Protein AAI columns: {len(prot_cols)}")

    # ------------------------------------------
    # 1. TAXONOMIC TREE AND TAXONOMIC DATASETS
    # ------------------------------------------

    tax_tree = None

    print("\n--- Generating Taxonomic Tree & Taxonomic Datasets ---")

    if not args.skip_taxonomy:
        tax_tree = build_taxonomy_tree(csv_path, tax_dir, DEFAULT_SLEEP)
    else:
        print("[INFO] Skipping taxonomy tree generation.")

    df_tax = map_taxonomic(df, best_col, tot_col, prot_cols)

    write_datasets(df_tax, tax_dir, best_col, tot_col, prot_cols)
    write_taxon_label_dataset(df_tax, tax_dir, include_taxid_in_label=args.label_taxid)

    print(f"[OK] Taxonomic datasets saved to: {tax_dir}")
    print(f"[OK] Taxonomic dataset rows/leaves: {len(df_tax)}")
    print("[INFO] Upload taxonomy datasets ONLY to taxonomy_tree_full.nwk.")
    print("[INFO] Upload itol_TAXON_LABELS.txt to display species names on taxid leaves.")

    # ------------------------------------------
    # 2. PHYLOGENETIC DATASETS
    # ------------------------------------------

    print("\n--- Generating Phylogenetic Datasets ---")

    df_phylo = map_phylogenetic(df, args.mode)
    write_datasets(df_phylo, phy_dir, best_col, tot_col, prot_cols)

    print(f"[OK] Phylogenetic datasets saved to: {phy_dir}")
    print(f"[OK] Phylogenetic dataset rows/leaves: {len(df_phylo)}")
    print("[INFO] Upload phylogenetic datasets ONLY to the matching Sourmash/local phylogeny tree.")

    # ------------------------------------------
    # 3. VALIDATION
    # ------------------------------------------

    if args.validate:
        if tax_tree is None:
            tax_tree = tax_dir / "taxonomy_tree_full.nwk"

        validate_dataset_ids_against_tree(tax_tree, tax_dir)

    print("\n=== Done ===")
    print("Important:")
    print("1. Use files from itol_out/taxonomic only with itol_out/taxonomic/taxonomy_tree_full.nwk.")
    print("2. Use files from itol_out/phylogenetic only with your matching phylogenetic tree.")
    print("3. Taxonomic IDs remain stable taxid-only IDs, e.g. taxid562.")
    print("4. Visible leaf names are controlled by itol_TAXON_LABELS.txt.")


if __name__ == "__main__":
    main()