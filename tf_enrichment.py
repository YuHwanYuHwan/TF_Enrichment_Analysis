#!/usr/bin/env python3
"""TF-target enrichment with an empirical null, tested separately for each reference DB.

Draw as many genes as the input from the background without replacement, count each TF's targets
among them, repeat n_perm times, p = #(null hits >= observed hits) / n_perm.
No +1, no multiple-testing correction.
See README.md for inputs and outputs.
"""
import argparse
import hashlib
import os
import platform
import shlex
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix

REF_COLS = ["TF_Ensembl", "TF_Entrez", "TF_symbol", "Target_Ensembl", "Target_Entrez", "Target_symbol"]
MAP_COLS = ["gene_id", "gene_symbol"]
EXCLUDED_COLS = ["set", "file", "id", "reason", "candidates"]
CHUNK_BYTES = 256 * 2**20  # caps the per-block arrays (draws, sparse sample matrix, null counts)


class AnalysisStop(Exception):
    """The data cannot be tested; the message is written to the log."""


def read_ids(path):
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    ids = (s.strip() for s in lines)
    return list(dict.fromkeys(s for s in ids if s and not s.startswith("#")))


def read_tsv(path, required):
    # NA markers ("NA" from R, "NULL", "#N/A", ...) are missing values, never gene symbols.
    # index_col=False: rows ending in an extra tab must not shift every column.
    df = pd.read_csv(path, sep="\t", dtype=str, index_col=False, encoding="utf-8-sig").fillna("")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    if df.apply(lambda col: col.str.contains("\n", regex=False)).to_numpy().any():
        raise ValueError(f"{path}: a quoted field spans several lines (unbalanced double quote?)")
    return df[required].apply(lambda col: col.str.strip())


def load_id_map(path):
    # Long format: one ID per row, and an ID may repeat; a repeated identical symbol counts once.
    m = read_tsv(path, MAP_COLS)
    known, symbols = set(), {}
    for gid, sym in zip(m["gene_id"], m["gene_symbol"]):
        if gid:
            known.add(gid)
            if sym:
                symbols.setdefault(gid, set()).add(sym)
    return known, symbols


def map_ids(ids, known, symbols, kind, source):
    mapped, excluded = set(), []
    for gid in ids:
        cands = symbols.get(gid, set())
        if gid not in known:
            reason = "ID_NOT_FOUND"
        elif not cands:
            reason = "SYMBOL_MISSING"
        elif len(cands) > 1:
            reason = "AMBIGUOUS_SYMBOL"
        else:
            mapped.add(next(iter(cands)))
            continue
        excluded.append({"set": kind, "file": str(source), "id": gid, "reason": reason,
                         "candidates": ";".join(sorted(cands))})
    return sorted(mapped), excluded


def load_reference(path):
    ref = read_tsv(path, REF_COLS)
    blank = (ref["TF_symbol"] == "") | (ref["Target_symbol"] == "")
    pairs = ref.loc[~blank, ["TF_symbol", "Target_symbol"]].drop_duplicates().reset_index(drop=True)
    return pairs, int(blank.sum()), len(ref)


def enrich_db(pairs, input_syms, bg_syms, n_perm, seed):
    """One DB. pairs: unique (TF_symbol, Target_symbol); input_syms must be a subset of bg_syms."""
    tfs = np.array(sorted(pairs["TF_symbol"].unique()), dtype=object)
    if not len(tfs):
        return pd.DataFrame({"TF": [], "p_value": [], "hit_genes": []})
    bg_pos = pd.Series(np.arange(len(bg_syms)), index=bg_syms)
    tf_pos = pd.Series(np.arange(len(tfs)), index=tfs)

    # Genes that are no TF's target add 0 to every count, so only target genes get a row.
    in_bg = pairs[pairs["Target_symbol"].isin(bg_pos.index)]
    gene = bg_pos.loc[in_bg["Target_symbol"]].to_numpy()
    tf = tf_pos.loc[in_bg["TF_symbol"]].to_numpy()
    target_genes = np.unique(gene)
    row_of = np.full(len(bg_syms), -1, dtype=np.int64)
    row_of[target_genes] = np.arange(len(target_genes))
    member = np.zeros((len(target_genes), len(tfs)), dtype=np.float32)
    member[row_of[gene], tf] = 1.0

    rows = row_of[bg_pos.loc[input_syms].to_numpy()]
    observed = member[rows[rows >= 0]].sum(axis=0)

    # Same seed for every DB: identical draws, so a DB's p-values do not depend on which other DBs run.
    rng = np.random.default_rng(seed)
    n_bg, n_in = len(bg_syms), len(input_syms)
    chunk = max(1, min(n_perm, CHUNK_BYTES // (8 * (len(tfs) + 5 * n_in))))
    exceed = np.zeros(len(tfs), dtype=np.int64)
    for start in range(0, n_perm, chunk):
        size = min(chunk, n_perm - start)
        draws = np.stack([rng.choice(n_bg, n_in, replace=False, shuffle=False) for _ in range(size)])
        hit_rows = row_of[draws]
        keep = hit_rows >= 0
        indptr = np.concatenate([[0], np.cumsum(keep.sum(axis=1))])
        sampled = csr_matrix((np.ones(indptr[-1], dtype=np.float32), hit_rows[keep], indptr),
                             shape=(size, len(target_genes)))
        exceed += (sampled @ member >= observed).sum(axis=0)

    hits = pairs[pairs["Target_symbol"].isin(set(input_syms))]
    hit_genes = hits.groupby("TF_symbol")["Target_symbol"].agg(lambda s: ",".join(sorted(s)))
    return pd.DataFrame({"TF": tfs, "p_value": exceed / n_perm,
                         "hit_genes": hit_genes.reindex(tfs).fillna("").to_numpy()})


def merge_results(results):
    merged = None
    for label, res in results.items():
        r = res.rename(columns={"p_value": f"{label}_p_value", "hit_genes": f"{label}_hit_genes"})
        merged = r if merged is None else merged.merge(r, on="TF", how="outer")
    return merged.sort_values("TF").reset_index(drop=True)


def run(input_path, background_path, id_map_path, ref_paths, n_perm, seed, log, excluded):
    """Returns ({db: per-DB table}, merged table). Appends to the caller's log/excluded lists."""
    if n_perm < 1 or seed < 0:
        raise ValueError("--n-perm must be >= 1 and --seed >= 0")
    labels = [Path(p).stem for p in ref_paths]
    if len({label.lower() for label in labels}) != len(labels):
        raise ValueError(f"reference file names must be unique, ignoring case (DB name = file name stem): {labels}")

    known, symbols = load_id_map(id_map_path)
    sets = {}
    for kind, path in (("input", input_path), ("background", background_path)):
        ids = read_ids(path)
        syms, exc = map_ids(ids, known, symbols, kind, path)
        excluded += exc
        sets[kind] = syms
        reasons = ", ".join(f"{r} {n}" for r, n in sorted(Counter(e["reason"] for e in exc).items())) or "none"
        log.append(f"{kind}: {len(ids)} unique IDs, excluded {len(exc)} ({reasons}), "
                   f"{len(ids) - len(exc) - len(syms)} collapsed into shared symbols, "
                   f"{len(syms)} unique symbols used ({path})")
    for kind, syms in sets.items():
        if not syms:
            raise AnalysisStop(f"{kind} is empty after ID conversion")
    outside = sorted(set(sets["input"]) - set(sets["background"]))
    if outside:
        raise AnalysisStop(f"{len(outside)} input symbols are not in the background, so the null could never "
                           f"draw them: {', '.join(outside)}")

    bg_set = set(sets["background"])
    results = {}
    for label, path in zip(labels, ref_paths):
        pairs, n_blank, n_rows = load_reference(path)
        targets = pairs["Target_symbol"].drop_duplicates()
        log.append(f"{label}: {n_rows} rows, {n_blank} with a blank symbol dropped, {len(pairs)} unique TF-target pairs, "
                   f"{pairs['TF_symbol'].nunique()} TFs, {len(targets)} targets "
                   f"({int((~targets.isin(bg_set)).sum())} not in background) ({path})")
        if pairs.empty:
            log.append(f"{label}: no usable TF-target pairs, DB skipped")
            continue
        results[label] = enrich_db(pairs, sets["input"], sets["background"], n_perm, seed)
    if not results:
        raise AnalysisStop("no reference DB had usable TF-target pairs")
    return results, merge_results(results)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same_file(a, b):
    try:
        return os.path.samefile(a, b)
    except OSError:  # a missing file is reported as not found further down
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description="TF-target enrichment with an empirical null, one test per TF per reference DB.")
    ap.add_argument("--input", required=True, help="gene IDs to test (Ensembl or Entrez), one per line")
    ap.add_argument("--background", required=True, help="background gene IDs, one per line")
    ap.add_argument("--ref-dir", default="reference",
                    help="folder of one species holding the id-map and one *.tsv per DB (default: %(default)s)")
    ap.add_argument("--id-map",
                    help="TSV with columns gene_id, gene_symbol (same species; other columns ignored) "
                         "(default: <ref-dir>/gene_id_symbol_master.tsv)")
    ap.add_argument("--ref", action="append",
                    help="TF-target reference TSV; repeat once per DB (DB name = file name stem). "
                         "Default: every *.tsv in --ref-dir except the id-map")
    ap.add_argument("--outdir", default="output", help="output folder (default: %(default)s)")
    ap.add_argument("--n-perm", type=int, default=100_000, help="number of random draws (default: %(default)s)")
    ap.add_argument("--seed", type=int, default=42, help="random seed, non-negative integer (default: %(default)s)")
    args = ap.parse_args(argv)
    if args.id_map is None:
        args.id_map = str(Path(args.ref_dir) / "gene_id_symbol_master.tsv")
    if args.ref:
        named = [r for r in args.ref if same_file(r, args.id_map)]
        if named:
            ap.error(f"--ref {', '.join(named)} is the id-map, not a TF-target reference")
    else:
        args.ref = sorted(str(p) for p in Path(args.ref_dir).glob("*.tsv") if not same_file(p, args.id_map))
        if not args.ref:
            ap.error(f"no reference files in {args.ref_dir}: put one *.tsv per DB there next to the id-map, "
                     "or pass --ref")
    # Reject bad arguments before --outdir is touched, so old results survive a typo.
    if args.n_perm < 1 or args.seed < 0:
        ap.error("--n-perm must be >= 1 and --seed >= 0")
    absent = [p for p in [args.input, args.background, args.id_map, *args.ref] if not Path(p).is_file()]
    if absent:
        ap.error(f"file not found: {', '.join(absent)}")

    out = Path(args.outdir)
    (out / "per_db").mkdir(parents=True, exist_ok=True)
    for old in [out / "merged.csv", *(out / "per_db").glob("*.csv")]:  # no stale results next to a new log
        old.unlink(missing_ok=True)
    log = [f"tf_enrichment {datetime.now().isoformat(timespec='seconds')}",
           f"command: {shlex.join(sys.argv if argv is None else ['tf_enrichment.py', *argv])}",
           f"python {platform.python_version()}, numpy {np.__version__}, pandas {pd.__version__}, scipy {scipy.__version__}",
           f"n_perm {args.n_perm}, seed {args.seed}"]
    excluded, status = [], 0
    try:
        for p in [args.input, args.background, args.id_map, *args.ref]:
            log.append(f"sha256 {sha256(p)} {p}")
        results, merged = run(args.input, args.background, args.id_map, args.ref,
                              args.n_perm, args.seed, log, excluded)
        for label, res in results.items():
            res.sort_values(["p_value", "TF"]).to_csv(out / "per_db" / f"{label}.csv", index=False)
        merged.to_csv(out / "merged.csv", index=False)
        log.append("done")
    except AnalysisStop as e:
        log.append(f"STOPPED: {e}")
        print(f"stopped: {e}", file=sys.stderr)
        status = 2
    except Exception as e:
        log.append(f"FAILED: {type(e).__name__}: {e}")
        raise
    finally:
        pd.DataFrame(excluded, columns=EXCLUDED_COLS).to_csv(out / "excluded_ids.tsv", sep="\t", index=False)
        (out / "run.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    return status


if __name__ == "__main__":
    sys.exit(main())
