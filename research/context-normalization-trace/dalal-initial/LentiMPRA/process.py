"""Deterministic build of the four ccRE lentiMPRA parquet files from raw/*.xlsx.

Grain: one row per *design* construct (S3 large-scale per cell line; S10 joint).
Every design row is preserved; measurement/QC/modeling exclusions become nulls,
never dropped rows. Source row order is preserved.

XLSX is read with a minimal stdlib zip/xml streaming reader (copied inline) because
the source workbooks are >100 MB and openpyxl read_only is ~5-10x slower here; the
reader was cross-validated against openpyxl to machine precision.
"""
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow  # noqa: F401  (parquet engine; imported to fail loudly if absent)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common import sha256, write_json

RAW = str(ROOT / "data/raw")
OUT = str(ROOT / "LentiMPRA")
FILES = {
    "S3": "2023-03-03628C-Table S3 - large-scale lib design (1).xlsx",
    "S4": "2023-03-03628C-Table S4 - largescaledata (1).xlsx",
    "S8": "2023-03-03628C-Table S8 - folds and performance.xlsx",
    "S10": "2023-03-03628C-Table S10 - joint lib design.xlsx",
    "S11": "2023-03-03628C-Table S11 - jointdata.xlsx",
    "S12": "2023-03-03628C-Table S12 - folds and performance.xlsx",
}
ADAPT5 = "AGGACCGGATCAACT"
ADAPT3 = "CATTGCGTGAACCGA"
CELLS = ["HepG2", "K562", "WTC11"]
SHUFFLED_CATEGORY = "negative control, shuffled"
EXPECTED_ACTIVE = {
    "HepG2": {"promoter": (11367, 20816), "potential enhancer": (50714, 118433)},
    "K562": {"promoter": (15362, 29376), "potential enhancer": (69820, 169260)},
    "WTC11": {"promoter": (5038, 9964), "potential enhancer": (11861, 45942)},
}
EXPECTED_ACTIVE_THRESHOLDS = {"HepG2": -0.2465, "K562": -0.1176, "WTC11": -0.3284}

# --------------------------------------------------------------------------- #
# Minimal streaming XLSX reader (stdlib only)
# --------------------------------------------------------------------------- #
_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_COL = re.compile(r"^([A-Z]+)")


def _col_to_idx(ref):
    m = _COL.match(ref)
    s = 0
    for ch in m.group(1):
        s = s * 26 + (ord(ch) - 64)
    return s - 1


def _local(tag):
    return tag.rsplit("}", 1)[-1]


class Book:
    """Opens one workbook, caches shared strings, exposes sheets by name."""

    def __init__(self, tag):
        self.zf = zipfile.ZipFile(os.path.join(RAW, FILES[tag]))
        wb = ET.fromstring(self.zf.read("xl/workbook.xml"))
        rels = ET.fromstring(self.zf.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {}
        for rel in rels:
            tgt = rel.get("Target")
            if not tgt.startswith("/"):
                tgt = "xl/" + tgt.lstrip("./")
            rid_to_target[rel.get("Id")] = tgt
        self.sheets = {sh.get("name"): rid_to_target[sh.get(_REL + "id")]
                       for sh in wb.iter(_NS + "sheet")}
        self.shared = self._load_shared()

    def _load_shared(self):
        out = []
        if "xl/sharedStrings.xml" not in self.zf.namelist():
            return out
        parts, in_si = [], False
        with self.zf.open("xl/sharedStrings.xml") as f:
            for ev, el in ET.iterparse(f, events=("start", "end")):
                ln = _local(el.tag)
                if ev == "start" and ln == "si":
                    in_si, parts = True, []
                elif ev == "end":
                    if ln == "t" and in_si:
                        parts.append(el.text or "")
                    elif ln == "si":
                        out.append("".join(parts))
                        in_si = False
                        el.clear()
        return out

    def rows(self, sheet_name):
        """Yield (row_index:int, list-of-cell-values). Numbers come back as raw
        strings; shared/inline strings as text; blanks as None."""
        path = self.sheets[sheet_name]
        with self.zf.open(path) as f:
            cells, ctype, cidx, vtext, itext = None, None, None, None, None
            for ev, el in ET.iterparse(f, events=("start", "end")):
                ln = _local(el.tag)
                if ev == "start":
                    if ln == "row":
                        cells = {}
                    elif ln == "c":
                        ctype = el.get("t")
                        ref = el.get("r")
                        cidx = _col_to_idx(ref) if ref else None
                        vtext, itext = None, []
                else:
                    if ln == "v":
                        vtext = el.text
                    elif ln == "t":
                        itext.append(el.text or "")
                    elif ln == "c":
                        if ctype == "s":
                            val = self.shared[int(vtext)] if vtext is not None else None
                        elif ctype in ("inlineStr", "str"):
                            val = "".join(itext) if ctype == "inlineStr" else vtext
                        elif ctype == "b":
                            val = vtext == "1"
                        else:
                            val = vtext  # number: keep raw string
                        if cidx is None:
                            cidx = (max(cells) + 1) if cells else 0
                        cells[cidx] = val
                        el.clear()
                    elif ln == "row":
                        if cells:
                            mx = max(cells)
                            yield int(el.get("r")), [cells.get(i) for i in range(mx + 1)]
                        else:
                            yield int(el.get("r")), []
                        el.clear()


_HEADER_TOKENS = {"name", "fold", "replicate"}


def read_table(book, sheet_name):
    """Return (header_list, iterator_of_data_rows). Header is the first row whose
    first cell is a known token; data rows follow it (blank rows skipped)."""
    it = book.rows(sheet_name)
    buffered, header, hrow = [], None, None
    for ridx, row in it:
        first = str(row[0]).strip().lower() if row and row[0] is not None else ""
        if first in _HEADER_TOKENS:
            header, hrow = row, ridx
            break
        buffered.append((ridx, row))
    if header is None:  # fallback: first non-empty row is the header
        for ridx, row in buffered:
            if any(c not in (None, "") for c in row):
                header, hrow = row, ridx
                break

    def data():
        for ridx, row in book.rows(sheet_name):
            if ridx <= hrow:
                continue
            if any(c is not None and c != "" for c in row):
                yield row
    return header, data()


def hmap(header):
    return {h: i for i, h in enumerate(header) if h not in (None, "")}


# --------------------------------------------------------------------------- #
# Parsing / normalization helpers
# --------------------------------------------------------------------------- #
def to_int(x):
    if x is None or x == "" or x == "NA":
        return None
    return int(float(x))


def to_float(x):
    if x is None or x == "" or x == "NA":
        return None
    return float(x)


def norm_chrom(x):
    if x is None or x == "" or x == "NA":
        return None
    x = str(x)
    return x if x.startswith("chr") else "chr" + x


def norm_strand(x):
    if x is None or x == "" or x == "NA":
        return None
    return str(x)


# Design sheets write the 2nd coordinate as _[<c>]; data sheets as __<c>_ .
# Same construct; this one-time fallback maps a data name to its design name.
_ALIAS = re.compile(r"__(.+?)_(?=$|_Reversed:|:)")


def alias_normalize(name):
    return _ALIAS.sub(lambda m: "_[" + m.group(1) + "]", name, count=1)


def parse_eo(name, cell):
    """Return (element_id, orientation) using the proven cell-specific grammar:
    WTC11-origin uses _F/_R; HepG2/K562/other use _Reversed:."""
    wtc11_origin = cell == "WTC11" or (cell == "joint" and name.startswith("WTC11_"))
    if wtc11_origin:
        if name.endswith("_R"):
            return name[:-2], "reverse"
        if name.endswith("_F"):
            return name[:-2], "forward"
        return name, "forward"
    if name.endswith("_Reversed:"):
        return name[: -len("_Reversed:")], "reverse"
    return name, "forward"


def canon(name, design_set):
    """Exact match first, then the one-time alias fallback."""
    if name in design_set:
        return name
    return alias_normalize(name)


def trim_insert(seq):
    assert seq is not None and len(seq) == 230, f"sequence not 230 nt: {seq!r}"
    assert seq[:15] == ADAPT5 and seq[-15:] == ADAPT3, f"bad adapters: {seq[:15]}..{seq[-15:]}"
    return seq[15:-15]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_design(book, sheet):
    """Ordered list of design row dicts + the set of design names."""
    header, data = read_table(book, sheet)
    hm = hmap(header)
    seq_col = next(i for h, i in hm.items() if h and ("230nt" in h or "sequence" in h.lower()))
    rows, names = [], set()
    for row in data:
        def g(col):
            i = hm.get(col)
            return row[i] if i is not None and i < len(row) else None
        name = g("name")
        rows.append({
            "name": name,
            "category": g("category"),
            "hg19_chrom": norm_chrom(g("chr.hg19")),
            "hg19_start": to_int(g("start.hg19")),
            "hg19_end": to_int(g("stop.hg19")),
            "hg19_strand": norm_strand(g("str.hg19")),
            "hg38_chrom": norm_chrom(g("chr.hg38")),
            "hg38_start": to_int(g("start.hg38")),
            "hg38_end": to_int(g("stop.hg38")),
            "hg38_strand": norm_strand(g("str.hg38")),
            "sequence": trim_insert(row[seq_col] if seq_col < len(row) else None),
        })
        names.add(name)
    return rows, names


def load_replicate_level(book, sheet, design_set):
    """dict (canon_name, replicate) -> (log2_full, barcodes). Asserts uniqueness
    and that every measurement name maps into the design set."""
    header, data = read_table(book, sheet)
    hm = hmap(header)
    ci_rep, ci_name = hm["replicate"], hm["name"]
    ci_log2 = hm["log2(rna/dna)"]
    ci_bc = hm["num_observed_barcodes"]
    out, dups, unmapped = {}, 0, 0
    for row in data:
        nm = row[ci_name]
        cn = canon(nm, design_set)
        if cn not in design_set:
            unmapped += 1
        key = (cn, int(row[ci_rep]))
        if key in out:
            dups += 1
        out[key] = (float(row[ci_log2]), int(row[ci_bc]))
    assert dups == 0, f"{sheet}: {dups} duplicate (name,replicate) keys"
    assert unmapped == 0, f"{sheet}: {unmapped} measurement names not mapped to design"
    return out


def load_summary(book, sheet, design_set):
    """dict canon_name -> {'reps':[r1,r2,r3], 'mean':m} (published, rounded)."""
    header, data = read_table(book, sheet)
    hm = hmap(header)
    ci = [hm[f"replicate {k} [log2(rna/dna)]"] for k in (1, 2, 3)]
    ci_mean = hm["mean"]
    out = {}
    for row in data:
        cn = canon(row[hm["name"]], design_set)
        assert cn not in out, f"{sheet}: duplicate summary name after normalization: {cn}"
        out[cn] = {"reps": [to_float(row[c]) if c < len(row) else None for c in ci],
                   "mean": to_float(row[ci_mean]) if ci_mean < len(row) else None}
    return out


def load_s8(book, cell):
    header, data = read_table(book, f"Fold assignment, {cell}")
    hm = hmap(header)
    ci_fold, ci_id, ci_obs = hm["Fold"], hm["Sequence ID"], hm["Observed log2(RNA/DNA)"]
    return {row[ci_id]: (float(row[ci_obs]), int(row[ci_fold])) for row in data}


def load_forward_reverse(book, cell):
    header, data = read_table(book, f"{cell}_forward_reverse")
    hm = hmap(header)
    ci_n, ci_f, ci_r = hm["name"], hm["forward [log2(rna/dna)]"], hm["reverse [log2(rna/dna)]"]
    return {row[ci_n]: (to_float(row[ci_f]), to_float(row[ci_r])) for row in data}


def load_benchmark_results(book, dataset):
    records = []
    cell = None
    models = None
    for _, row in book.rows("Benchmarking results"):
        first = row[0] if row else None
        if first in CELLS:
            cell = first
            models = None
        elif str(first).strip().lower() == "test fold":
            models = [str(value) for value in row[1:] if value not in (None, "")]
        elif cell and models and first not in (None, ""):
            fold = to_int(first)
            if fold is None:
                continue
            for model, value in zip(models, row[1:]):
                records.append(
                    {
                        "dataset": dataset,
                        "cell_line": cell,
                        "test_fold": fold,
                        "model": model,
                        "pearson_r": float(value),
                    }
                )
    expected_models = 6 if dataset == "large_scale" else 5
    check(
        len(records) == len(CELLS) * 10 * expected_models,
        f"{dataset}: unexpected published benchmark row count",
    )
    return pd.DataFrame(
        {
            "dataset": pd.array([row["dataset"] for row in records], dtype="string"),
            "cell_line": pd.array([row["cell_line"] for row in records], dtype="string"),
            "test_fold": pd.array([row["test_fold"] for row in records], dtype="Int8"),
            "model": pd.array([row["model"] for row in records], dtype="string"),
            "pearson_r": pd.Series(
                [row["pearson_r"] for row in records],
                dtype="float64",
            ),
        }
    )


# --------------------------------------------------------------------------- #
# Frame builders
# --------------------------------------------------------------------------- #
IDENTITY_COLS = ["name", "element_id", "category", "sequence"]
COORDINATE_COLS = [
    "hg19_chrom", "hg19_start", "hg19_end", "hg19_strand",
    "hg38_chrom", "hg38_start", "hg38_end", "hg38_strand",
    "orientation",
]
LARGE_COLS = [
    *IDENTITY_COLS,
    *COORDINATE_COLS,
    "replicate_1_log2_rna_dna", "replicate_1_num_observed_barcodes",
    "replicate_2_log2_rna_dna", "replicate_2_num_observed_barcodes",
    "replicate_3_log2_rna_dna", "replicate_3_num_observed_barcodes",
    "mean_log2_rna_dna", "observed_log2_rna_dna",
    "active", "fold",
]
STR_COLS = {"name", "element_id", "category", "orientation",
            "hg19_chrom", "hg19_strand", "hg38_chrom", "hg38_strand", "sequence"}
INT_COLS = {"hg19_start", "hg19_end", "hg38_start", "hg38_end", "fold"}


def frame_from_records(records, columns):
    """Build a DataFrame with explicit nullable dtypes and exact column order."""
    data = {}
    for c in columns:
        vals = [r.get(c) for r in records]
        if c in STR_COLS:
            data[c] = pd.array(vals, dtype="string")
        elif c == "active" or c.endswith("_active"):
            data[c] = pd.array(vals, dtype="Int8")
        elif c in INT_COLS or c.endswith("num_observed_barcodes"):
            data[c] = pd.array([None if v is None else int(v) for v in vals], dtype="Int64")
        else:
            data[c] = pd.Series([np.nan if v is None else float(v) for v in vals], dtype="float64")
    return pd.DataFrame(data, columns=columns)


def add_active_labels(records, mean_col, active_col):
    """Call measured rows active above the shuffled-control 95th percentile."""
    shuffled = [
        row[mean_col]
        for row in records
        if row["category"] == SHUFFLED_CATEGORY and row.get(mean_col) is not None
    ]
    check(shuffled, f"{active_col}: no measured shuffled controls")
    threshold = float(pd.Series(shuffled, dtype="float64").quantile(0.95, interpolation="linear"))
    for row in records:
        mean = row.get(mean_col)
        row[active_col] = None if mean is None else int(mean >= threshold)
    return threshold, len(shuffled)


def build_large_scale(cell, books):
    design, design_set = load_design(books["S3"], f"{cell} large-scale")
    rep = load_replicate_level(books["S4"], f"{cell}_replicate_level", design_set)
    summ = load_summary(books["S4"], f"{cell}_summary_data", design_set)
    s8 = load_s8(books["S8"], cell)
    records = []
    for d in design:
        eid, orient = parse_eo(d["name"], cell)
        rec = dict(d)
        rec["element_id"], rec["orientation"] = eid, orient
        for k in (1, 2, 3):
            lb = rep.get((d["name"], k))
            rec[f"replicate_{k}_num_observed_barcodes"] = None if lb is None else lb[1]
            rec[f"replicate_{k}_log2_rna_dna"] = None if lb is None else lb[0]
        rec["mean_log2_rna_dna"] = (summ.get(d["name"]) or {}).get("mean")
        obs = s8.get(eid)  # S8 is element-level; duplicated onto both orientations
        rec["observed_log2_rna_dna"] = None if obs is None else obs[0]
        rec["fold"] = None if obs is None else obs[1]
        records.append(rec)
    threshold, n_shuffled = add_active_labels(records, "mean_log2_rna_dna", "active")
    return frame_from_records(records, LARGE_COLS), design_set, rep, summ, s8, threshold, n_shuffled


def joint_columns():
    cols = [*IDENTITY_COLS, *COORDINATE_COLS]
    for p in ("hepg2", "k562", "wtc11"):
        for k in (1, 2, 3):
            cols += [
                f"{p}_replicate_{k}_log2_rna_dna",
                f"{p}_replicate_{k}_num_observed_barcodes",
            ]
        cols += [
            f"{p}_mean_log2_rna_dna",
            f"{p}_element_specificity_score",
            f"{p}_observed_element_specificity_score",
        ]
    cols += [f"{p}_active" for p in ("hepg2", "k562", "wtc11")]
    cols += ["fold"]
    return cols


def build_joint(books):
    design, design_set = load_design(books["S10"], "joint library")
    s11 = books["S11"]
    rep = {c: load_replicate_level(s11, f"{c}_replicate_level", design_set) for c in CELLS}
    summ = {c: load_summary(s11, f"{c}_summary_data", design_set) for c in CELLS}
    # published ESS (all_cell_types_ESS) and observed ESS / fold (S12)
    h_ess, d_ess = read_table(s11, "all_cell_types_ESS")
    hm = hmap(h_ess)
    ci_ess = {c: hm[f"ElementSpecificityScore_{c}"] for c in CELLS}
    ess = {}
    for row in d_ess:
        ess[canon(row[hm["name"]], design_set)] = {c: to_float(row[ci_ess[c]]) for c in CELLS}
    h12, d12 = read_table(books["S12"], "Fold assignment, joint library")
    hm12 = hmap(h12)
    ci12 = {c: hm12[f"{c}, Element Specificity Score"] for c in CELLS}
    s12 = {}
    for row in d12:
        cn = canon(row[hm12["Sequence ID"]], design_set)
        s12[cn] = ({c: to_float(row[ci12[c]]) for c in CELLS}, int(row[hm12["Fold"]]))

    cols = joint_columns()
    records = []
    for d in design:
        eid, orient = parse_eo(d["name"], "joint")
        rec = dict(d)
        rec["element_id"], rec["orientation"] = eid, orient
        for p, c in zip(("hepg2", "k562", "wtc11"), CELLS):
            for k in (1, 2, 3):
                lb = rep[c].get((d["name"], k))
                rec[f"{p}_replicate_{k}_num_observed_barcodes"] = None if lb is None else lb[1]
                rec[f"{p}_replicate_{k}_log2_rna_dna"] = None if lb is None else lb[0]
            rec[f"{p}_mean_log2_rna_dna"] = (summ[c].get(d["name"]) or {}).get("mean")
            rec[f"{p}_element_specificity_score"] = (ess.get(d["name"]) or {}).get(c)
            obs = s12.get(d["name"])
            rec[f"{p}_observed_element_specificity_score"] = None if obs is None else obs[0][c]
        rec["fold"] = None if d["name"] not in s12 else s12[d["name"]][1]
        records.append(rec)

    thresholds = {}
    for p in ("hepg2", "k562", "wtc11"):
        thresholds[p] = add_active_labels(
            records,
            mean_col=f"{p}_mean_log2_rna_dna",
            active_col=f"{p}_active",
        )
    frame = frame_joint(records, cols)
    return frame, design_set, rep, summ, ess, s12, thresholds


def frame_joint(records, columns):
    data = {}
    for c in columns:
        vals = [r.get(c) for r in records]
        if c in STR_COLS:
            data[c] = pd.array(vals, dtype="string")
        elif c == "active" or c.endswith("_active"):
            data[c] = pd.array(vals, dtype="Int8")
        elif c == "fold" or c.endswith("num_observed_barcodes"):
            data[c] = pd.array([None if v is None else int(v) for v in vals], dtype="Int64")
        elif c in ("hg19_start", "hg19_end", "hg38_start", "hg38_end"):
            data[c] = pd.array([None if v is None else int(v) for v in vals], dtype="Int64")
        else:
            data[c] = pd.Series([np.nan if v is None else float(v) for v in vals], dtype="float64")
    return pd.DataFrame(data, columns=columns)


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
def check(cond, msg):
    if not cond:
        raise AssertionError("VALIDATION FAILED: " + msg)


def validate_common(df, n_expected, design_order, label):
    check(len(df) == n_expected, f"{label}: row count {len(df)} != {n_expected}")
    check(df["name"].is_unique, f"{label}: names not unique")
    check(list(df["name"]) == design_order, f"{label}: row order not preserved")
    seqs = df["sequence"].dropna()
    check((seqs.str.len() == 200).all(), f"{label}: sequence not 200 nt")
    alpha = set("".join(seqs.tolist()))
    # designed constructs may contain IUPAC ambiguity codes, not only A/C/G/T
    check(alpha <= set("ACGTNRYSWKMBDHV"), f"{label}: non-IUPAC sequence alphabet {alpha}")
    folds = df["fold"].dropna().unique().tolist()
    check(set(folds) <= set(range(1, 11)), f"{label}: folds out of 1..10: {sorted(folds)}")


def validate_summary_rules(rep, summ, cell, design_set):
    """Summary rep columns are round(full-precision replicate log2, 3) and are
    non-null exactly for the kept replicates (bc>=10 for S4; all present for S11).
    The published mean is round(mean(FULL-precision kept reps), 3) -- the mean uses
    full precision even though the individual rep columns are rounded.
    Summary-only aggregate rows (e.g. 'no_BC') that are not design elements and
    carry no replicate_level data are not part of the output and are skipped."""
    for name, s in summ.items():
        if name not in design_set:
            continue
        kept_full = []
        for k in (1, 2, 3):
            v = s["reps"][k - 1]
            raw = rep.get((name, k))
            if v is not None:
                check(raw is not None, f"{cell}: summary rep present but no replicate_level row {name} r{k}")
                check(abs(v - round(raw[0], 3)) < 1e-9,
                      f"{cell}: summary rep {name} r{k} {v} != round(raw,3) {round(raw[0],3)}")
                kept_full.append(raw[0])
        if s["mean"] is not None and kept_full:
            check(abs(s["mean"] - round(float(np.mean(kept_full)), 3)) < 1e-9,
                  f"{cell}: mean {s['mean']} != round(mean(full kept),3) for {name}")
        check((s["mean"] is not None) == bool(kept_full),
              f"{cell}: summary mean presence does not match retained replicates for {name}")


def validate_s4_bc_filter(rep, summ, cell):
    """S4 large-scale: summary rep non-null IFF barcode>=10."""
    bad = 0
    for (name, k), (log2, bc) in rep.items():
        present = name in summ and summ[name]["reps"][k - 1] is not None
        if (bc >= 10) != present:
            bad += 1
    check(bad == 0, f"{cell}: {bad} replicate rows violate bc>=10 <-> summary-present")


def validate_joint_replicates(rep, summ, cell, design_set):
    """Joint summaries retain every replicate-level measurement."""
    missing = 0
    for name, k in rep:
        if name in design_set and (name not in summ or summ[name]["reps"][k - 1] is None):
            missing += 1
    check(missing == 0, f"joint:{cell}: {missing} replicate rows absent from summary")


def validate_active_labels(df, mean_col, active_col, threshold, label):
    measured = df[mean_col].notna()
    check(df.loc[~measured, active_col].isna().all(),
          f"{label}: unmeasured rows must have null {active_col}")
    expected = (df.loc[measured, mean_col] >= threshold).astype("int8").to_numpy()
    observed = df.loc[measured, active_col].astype("int8").to_numpy()
    check(np.array_equal(observed, expected), f"{label}: active labels do not match threshold")
    check(set(df[active_col].dropna().unique()) <= {0, 1},
          f"{label}: {active_col} contains values outside 0/1")


def validate_paper_active_counts(df, cell, threshold, n_shuffled):
    check(abs(threshold - EXPECTED_ACTIVE_THRESHOLDS[cell]) < 1e-12,
          f"{cell}: active threshold {threshold} != {EXPECTED_ACTIVE_THRESHOLDS[cell]}")
    for category, (expected_active, expected_total) in EXPECTED_ACTIVE[cell].items():
        rows = df[(df["category"] == category) & df["mean_log2_rna_dna"].notna()]
        check(len(rows) == expected_total,
              f"{cell} {category}: denominator {len(rows)} != {expected_total}")
        check(int(rows["active"].sum()) == expected_active,
              f"{cell} {category}: active {int(rows['active'].sum())} != {expected_active}")
    check(n_shuffled == {"HepG2": 311, "K562": 247, "WTC11": 83}[cell],
          f"{cell}: measured shuffled controls {n_shuffled} unexpected")


def validate_s8(books, dfs):
    maxerr = 0.0
    for cell in CELLS:
        fr = load_forward_reverse(books["S4"], cell)
        s8 = load_s8(books["S8"], cell)
        act = {e: (fr[e][0] + fr[e][1]) / 2 for e in s8 if e in fr and None not in fr[e]}
        check(len(act) == len(s8), f"S8 {cell}: {len(s8)-len(act)} elements missing forward/reverse")
        v = np.array(list(act.values()))
        mu, sd = v.mean(), v.std(ddof=1)
        err = max(abs((act[e] - mu) / sd - s8[e][0]) for e in act)
        maxerr = max(maxerr, err)
        check(err <= 1e-12, f"S8 {cell}: zscore reproduction err {err:.2e} > 1e-12")
    return maxerr


def validate_joint_ess(books, ess, s12, design_set):
    s11 = books["S11"]
    # all_cell_types_summary agrees with per-cell summary means
    h, d = read_table(s11, "all_cell_types_summary")
    hm = hmap(h)
    ci = {c: hm[f"{c} [log2(rna/dna)]"] for c in CELLS}
    acts = {}
    for row in d:
        acts[canon(row[hm["name"]], design_set)] = {c: to_float(row[ci[c]]) for c in CELLS}
    # ESS sums to zero and multiples of 1/18000
    esum, e18 = 0.0, 0.0
    for name, e in ess.items():
        vals = [e[c] for c in CELLS]
        if all(v is not None for v in vals):
            esum = max(esum, abs(sum(vals)))
            for v in vals:
                e18 = max(e18, abs(v * 18000 - round(v * 18000)))
    check(esum < 1e-9, f"S11 ESS sum-to-zero violated: {esum:.2e}")
    check(e18 < 1e-9, f"S11 ESS not multiples of 1/18000: {e18:.2e}")
    # S12 observed ESS == per-cell zscore(ddof=1) of S11 ESS over the S12 subset
    s12_names = list(s12)
    maxerr = 0.0
    for c in CELLS:
        base = np.array([ess[n][c] for n in s12_names])
        mu, sd = base.mean(), base.std(ddof=1)
        err = max(abs((ess[n][c] - mu) / sd - s12[n][0][c]) for n in s12_names)
        maxerr = max(maxerr, err)
        check(err <= 1e-12, f"S12 {c}: zscore reproduction err {err:.2e} > 1e-12")
    return acts, esum, e18, maxerr


def validate_readback(df, path, label):
    rb = pd.read_parquet(path, engine="pyarrow")
    pd.testing.assert_frame_equal(df, rb, check_dtype=True)
    print(f"    readback equal (dtypes/nulls): {label}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def summarize(df, label):
    nn = {c: int(df[c].notna().sum()) for c in df.columns}
    key = ["name", "active", "observed_log2_rna_dna", "fold"] if "active" in df.columns \
        else ["name", "hepg2_active", "hepg2_observed_element_specificity_score", "fold"]
    print(f"  [{label}] rows={len(df)} cols={len(df.columns)} "
          f"non-null(sample): " + ", ".join(f"{k}={nn[k]}" for k in key))


def main():
    os.makedirs(OUT, exist_ok=True)
    books = {t: Book(t) for t in FILES}
    print("Loaded workbooks:", ", ".join(FILES))

    # ---- large-scale, per cell line ----
    n_expected = {"HepG2": 164307, "K562": 243780, "WTC11": 75542}
    outputs = {}
    for cell in CELLS:
        print(f"\n== building large-scale {cell} ==")
        df, design_set, rep, summ, s8, active_threshold, n_shuffled = build_large_scale(cell, books)
        design_order = list(df["name"])
        check(len(design_set) == len(design_order), f"{cell}: design names not unique")
        check(all(alias_normalize(n) == n for n in design_set),
              f"{cell}: alias fallback alters an already-matched design name")
        validate_common(df, n_expected[cell], design_order, cell)
        validate_summary_rules(rep, summ, cell, design_set)
        validate_s4_bc_filter(rep, summ, cell)
        validate_active_labels(
            df,
            mean_col="mean_log2_rna_dna",
            active_col="active",
            threshold=active_threshold,
            label=cell,
        )
        validate_paper_active_counts(df, cell, active_threshold, n_shuffled)
        print(f"  active threshold={active_threshold:.6f} from n={n_shuffled} shuffled controls")
        outputs[cell] = df
        summarize(df, cell)

    s8_err = validate_s8(books, outputs)
    print(f"\nS8 target reproduction max error (all cells): {s8_err:.2e}")

    # ---- joint ----
    print("\n== building joint ==")
    jdf, jset, jrep, jsumm, jess, js12, joint_thresholds = build_joint(books)
    jorder = list(jdf["name"])
    check(len(jset) == len(jorder), "joint: design names not unique")
    check(all(alias_normalize(n) == n for n in jset),
          "joint: alias fallback alters an already-matched design name")
    validate_common(jdf, 60000, jorder, "joint")
    for c in CELLS:
        validate_summary_rules(jrep[c], jsumm[c], f"joint:{c}", jset)
        validate_joint_replicates(jrep[c], jsumm[c], c, jset)
    for p in ("hepg2", "k562", "wtc11"):
        threshold, n_shuffled = joint_thresholds[p]
        validate_active_labels(
            jdf,
            mean_col=f"{p}_mean_log2_rna_dna",
            active_col=f"{p}_active",
            threshold=threshold,
            label=f"joint:{p}",
        )
        print(f"  joint {p} active threshold={threshold:.6f} "
              f"from n={n_shuffled} shuffled controls")
    acts, esum, e18, s12_err = validate_joint_ess(books, jess, js12, jset)
    # per-cell summary mean agrees with all_cell_types_summary
    agree = 0
    for c in CELLS:
        for name, m in ((n, s["mean"]) for n, s in jsumm[c].items()):
            if name not in jset:
                continue
            a = acts.get(name, {}).get(c)
            if m is not None and a is not None:
                check(abs(m - a) < 1e-9, f"joint {c}: per-cell mean {m} != all_cell_types {a} for {name}")
                agree += 1
    print(f"joint: ESS sum|max={esum:.2e}, 1/18000 residual={e18:.2e}, "
          f"all_cell_types==per_cell means checked={agree}")
    print(f"S12 observed-ESS reproduction max error (all cells): {s12_err:.2e}")
    outputs["joint"] = jdf
    summarize(jdf, "joint")
    outputs["large_scale_performance"] = load_benchmark_results(
        books["S8"],
        dataset="large_scale",
    )
    outputs["joint_performance"] = load_benchmark_results(
        books["S12"],
        dataset="joint",
    )

    # ---- write + readback ----
    print("\n== writing parquet + readback validation ==")
    paths = {
        "HepG2": "lentimpra-hepg2.parquet",
        "K562": "lentimpra-k562.parquet",
        "WTC11": "lentimpra-wtc11.parquet",
        "joint": "lentimpra-joint.parquet",
        "large_scale_performance": "published-performance-large-scale.parquet",
        "joint_performance": "published-performance-joint.parquet",
    }
    for key, fn in paths.items():
        p = os.path.join(OUT, fn)
        temporary = p + ".tmp"
        outputs[key].to_parquet(temporary, engine="pyarrow", index=False)
        os.replace(temporary, p)
        validate_readback(outputs[key], p, fn)
        print(f"    wrote {fn}: {os.path.getsize(p)/1e6:.2f} MB, rows={len(outputs[key])}, "
              f"cols={len(outputs[key].columns)}")

    input_manifest = ROOT / "data/input_manifest.json"
    check(input_manifest.is_file(), "missing data/input_manifest.json")
    write_json(
        ROOT / "LentiMPRA/manifest.json",
        {
            "paper_doi": "10.1038/s41586-024-08430-9",
            "inputs": {
                "data/input_manifest.json": sha256(input_manifest),
                **{
                    f"data/raw/{filename}": sha256(Path(RAW) / filename)
                    for filename in FILES.values()
                },
            },
            "outputs": {
                filename: {
                    "rows": len(outputs[key]),
                    "columns": len(outputs[key].columns),
                    "sha256": sha256(ROOT / "LentiMPRA" / filename),
                }
                for key, filename in paths.items()
            },
        },
    )
    print("\nALL VALIDATIONS PASSED")


if __name__ == "__main__":
    main()
