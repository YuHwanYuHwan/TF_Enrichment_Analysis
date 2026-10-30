# TF Enrichment Analysis

Tests whether the targets of each transcription factor (TF) are over-represented in a gene list,
using an empirical (permutation) null. Each TF–target reference database (DB) is tested on its
own, and the per-DB results are then placed side by side in one table.

The tool is a single Python script that reads only local files: it does not download databases or
call ID-conversion services. Preparing the reference files is up to you (see [Inputs](#inputs)).

## Requirements

Python 3 with numpy, pandas and scipy. Tested on Python 3.12 (numpy 1.26, pandas 2.2, scipy 1.14)
and Python 3.13 (numpy 2.4, pandas 3.0, scipy 1.17). There is nothing to install beyond these
packages.

## Quick start

Set up a pipeline folder with one reference folder per species:

```
TF_Enrichment_Pipeline/
├── input/         gene list to test, background list
├── reference/
│   ├── mouse/     gene_id_symbol_master.tsv (the id-map) and one *.tsv per TF–target DB
│   └── human/     the same files for human
├── output/        results (created if missing)
└── script/        tf_enrichment.py
```

Then run from the pipeline folder, naming the species folder with `--ref-dir`:

```bash
cd TF_Enrichment_Pipeline
python script/tf_enrichment.py --ref-dir reference/mouse --input input/genes.txt --background input/background.txt
```

In the `--ref-dir` folder, `gene_id_symbol_master.tsv` is the id-map that converts IDs to symbols
and every other `*.tsv` file is tested as one DB; results go to `output/`. These defaults are
relative to the folder you run from; any of them can be replaced with an explicit path:

```bash
python tf_enrichment.py \
  --input genes.txt --background background.txt \
  --id-map id_map.tsv \
  --ref refs/DB_A.tsv --ref refs/DB_B.tsv \
  --outdir results/
```

When `--ref` is omitted, only the id-map in use is left out of the DBs. So if `--id-map` names
another file, a `gene_id_symbol_master.tsv` still in `--ref-dir` is read as a DB and the run fails
with a missing-column error; name the DBs with `--ref` as well, as in the example above.

| Option | Default | Meaning |
|---|---|---|
| `--input` | required | Gene IDs to test |
| `--background` | required | Gene IDs the null distribution draws from |
| `--ref-dir` | `reference` | Folder of one species holding the id-map and one `*.tsv` per DB |
| `--id-map` | `<ref-dir>/gene_id_symbol_master.tsv` | ID-to-symbol table for the same species |
| `--ref` | every `*.tsv` in `<ref-dir>` except the id-map | One TF–target reference file per DB; repeat the option for several DBs |
| `--outdir` | `output` | Output folder |
| `--n-perm` | 100000 | Number of random draws |
| `--seed` | 42 | Random seed (non-negative integer) |

Keep one species per reference folder (`--ref-dir`): all DBs in it are tested against the same
id-map and background, so a reference from another species would give p = 1 for all of its TFs.
Give each species its own `--outdir` as well (for example `output/mouse` and `output/human`),
because a run replaces the results already in its `--outdir`.

## Inputs

All files are UTF-8 text (a BOM is fine). Tables are tab-separated with a header row; input and
background are plain text files with one ID per line.

**`--input` and `--background`** — Ensembl gene IDs or Entrez Gene IDs (both may be mixed),
written exactly as in the id-map: a version suffix such as `ENSMUSG00000000001.4` does not match
`ENSMUSG00000000001`. Blank lines and lines starting with `#` are ignored, and a repeated ID counts
once. Every gene in the background is one the null may draw, and all DBs share the same
background, so choose it deliberately (for example the species' whole gene list, or only the genes
detected in your assay) and report the choice. The background must contain every input gene and
is normally much larger than the input; if both are the same list, every p is 1. The IDs below
show the format only:

```
ENSMUSG00000000001
ENSMUSG00000000028
67608
```

**`--id-map`** — a long table with one ID per row. Columns `gene_id` and `gene_symbol` are
required; other columns (such as `id_type` below) are ignored. Cover every gene of the species
(not only DB targets), with Ensembl and Entrez IDs on separate rows. One ID may appear on several
rows: rows repeating the same symbol count once, and an ID with two different symbols is excluded
as `AMBIGUOUS_SYMBOL`. The rows below show the format only:

```
id_type	gene_id	gene_symbol
ensembl	ENSMUSG00000000001	Gnai3
entrez	14679	Gnai3
```

**`--ref`** — one file per DB with the six columns below. One row is one TF–target relation. The
test uses only the two symbol columns: the four ID columns must be present but may be blank, and
they are not copied to the outputs. Rows with a blank TF or target symbol are dropped; a file left
with no pairs is skipped. The DB name used in the output is the file name without its extension
(`reference/mouse/DB_A.tsv` → `DB_A`). Only `.tsv` files other than the id-map are read from
`--ref-dir`, so other files there (notes, metadata) are ignored. The rows below show the format
only; they are not real regulatory relations.

```
TF_Ensembl	TF_Entrez	TF_symbol	Target_Ensembl	Target_Entrez	Target_symbol
ENSMUSG00000000078	23849	Klf6	ENSMUSG00000000001	14679	Gnai3
ENSMUSG00000000078	23849	Klf6	ENSMUSG00000000028	12544	Cdc45
```

### How IDs become genes

Input and background IDs are converted to symbols with the id-map, and every count is made on
unique symbols. Several IDs that map to the same symbol count as one gene. Matching is exact: the
script does not guess aliases or change letter case, so the reference symbols must already be the
official symbols of the same species and the same annotation release as the id-map.

**A mismatch does not stop the run.** A target whose symbol is not in the background can never be
hit, which pushes its TF's p towards 1; a human DB used with a mouse id-map gives p = 1 for every
TF. `run.log` reports, for each DB, how many targets are not in the background — check that count
after every run.

An ID that cannot be converted is dropped, the run continues, and the ID is listed in
`excluded_ids.tsv` with one of these reasons:

| Reason | Meaning |
|---|---|
| `ID_NOT_FOUND` | The ID is not in the id-map |
| `SYMBOL_MISSING` | The ID is in the id-map but has no symbol |
| `AMBIGUOUS_SYMBOL` | The ID maps to several different symbols; none is picked, and the candidates are listed |

Standard missing-value markers (`NA`, `NULL`, `#N/A`, …) in any table are read as blank, never as
gene symbols.

## How the test works

For each DB separately:

1. All rows of the same TF are merged into one target set; a repeated TF–target pair counts once.
2. **Observed** = number of the TF's targets among the input genes.
3. **Null**: draw as many genes as the input from the background, uniformly and without
   replacement, and count each TF's targets among them. One draw is shared by all TFs of the DB.
4. Repeat step 3 `--n-perm` times.
5. **p = (number of draws with hits ≥ observed) / n_perm**

Notes:

- There is no +1 correction, so `p = 0` means no draw reached the observed count; the smallest
  non-zero p is 1 / n_perm.
- No multiple-testing correction is applied.
- Every DB starts from the same seed, so the p-values of one DB do not change when other DBs are
  added to or removed from the run.
- A significant p means the TF's targets are over-represented in the input list. It does not say
  whether the TF is activated or repressed.

## Outputs

| File | Content |
|---|---|
| `per_db/<DB>.csv` | `TF`, `p_value`, `hit_genes` (input genes that are targets of the TF, comma-separated; `N/A` when none is), sorted by p-value |
| `merged.csv` | `TF`, then `<DB>_p_value` and `<DB>_hit_genes` for each DB in `--ref` order (a skipped DB has no columns). `N/A` in both columns means the TF is not in that DB; a p-value with `N/A` genes means no input gene is a target |
| `excluded_ids.tsv` | IDs dropped during conversion: input or background, source file, ID, reason, candidate symbols |
| `run.log` | Date, command, package versions, seed, n_perm, SHA-256 of every input file, gene counts after conversion, per-DB counts (rows, unique pairs, TFs, targets, targets outside the background) and the final status |

`merged.csv` holds each DB's own results side by side; no combined p-value, DB count or union of
hit genes is computed. Filter its p-value columns to pick TFs supported by the DBs you choose.

Reusing an `--outdir` first deletes its old `merged.csv` and `per_db/*.csv`, so every result file
belongs to the run that `run.log` describes.

## When a run stops

| Situation | Exit code | `run.log` ends with |
|---|---|---|
| Finished | 0 | `done` |
| Input or background empty after ID conversion | 2 | `STOPPED: ...` |
| An input gene is not in the background (the null could never draw it) | 2 | `STOPPED: ...` with the genes listed |
| No reference file has usable TF–target pairs | 2 | `STOPPED: ...` |
| Malformed file (missing column, unreadable encoding, a quoted field spanning several lines, duplicate DB names) | 1 | `FAILED: ...` |

Command-line errors (unknown option, missing required option, a value that is not an integer,
`--n-perm` below 1, a negative `--seed`, an input file that does not exist, the id-map given as a
`--ref`, no DB file in `--ref-dir` when `--ref` is omitted) exit with 2 before anything in
`--outdir` is created or deleted. A reference file without usable pairs is skipped and
noted in the log while the other DBs still run.

## Runtime

Measured on one CPU core with the default 100,000 draws and a 362-gene input list:

| Reference DB | Background | Time | Memory peak |
|---|---|---|---|
| Curated, 1,165 TFs, 43k pairs | 25k genes | about 8 s | about 0.5 GB |
| ChIP-scale, 1,500 TFs, 2.9M pairs | 55k genes | about 22 s | about 1.2 GB |

Time grows roughly in proportion to the input size and to `--n-perm`: the ChIP-scale DB with a
3,000-gene input took about 140 s.

## License

MIT — see [LICENSE](LICENSE).

## Author

YuHwan (yh9472@ajou.ac.kr)
