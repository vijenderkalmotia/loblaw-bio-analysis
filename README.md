# Immune cell-count analysis

`cell-count.csv` has 10,500 samples from 3,500 subjects across 3 projects. The CSV uses `sample` (not sample_id), `condition` (not indication), and `sex` (not gender). Cell counts are `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, and `monocyte`.

If another treatment such as quintazide showed up in the `treatment` column, it would load with the same schema. This file does not include quintazide.

## Run it (including GitHub Codespaces)

1. Open the repo in Codespaces (or clone it locally).
2. In the terminal:

```bash
make setup
make pipeline
make dashboard
```

`make setup` installs pandas, scipy, statsmodels, matplotlib, and Streamlit.

`make pipeline` rebuilds `cell_count.db` from the CSV, then writes every file under `outputs/`.

`make dashboard` starts Streamlit. In Codespaces, use **Ports** and open the forwarded Streamlit port (usually 8501). Locally it is [http://localhost:8501](http://localhost:8501).

If `python` is not on your PATH, the Makefile uses `python3`. Override with `make pipeline PYTHON=python`.

Rerunning `make pipeline` is safe: `load_data.py` drops and recreates the tables, so you do not get duplicate rows.

## Database schema

Four tables, matching how the assay is organized:

```text
projects (1) ──< subjects (1) ──< samples (1) ──< cell_counts
```

| Table | Grain | Columns |
| --- | --- | --- |
| `projects` | one study | `project_id` |
| `subjects` | one person | `subject_id`, `project_id`, `condition`, `age`, `sex`, `treatment`, `response` |
| `samples` | one blood draw | `sample_id`, `subject_id`, `sample_type`, `time_from_treatment_start` |
| `cell_counts` | one population in one sample | `sample_id`, `population`, `count` |

`response` is NULL for healthy donors (blank in the CSV). Foreign keys are on. Filter columns are indexed.

### Why this layout

In the CSV, each subject sits in exactly one project, and `condition`, `age`, `sex`, `treatment`, and `response` do not change across that subject's samples. Those fields belong on `subjects`, not copied onto every sample row.

Each sample is a visit (`PBMC` or `WB`, day 0 / 7 / 14). Counts are stored long, not as five extra columns on `samples`, so a new population is a new row instead of an `ALTER TABLE`. Frequency math is then `SUM(count)` per sample.

`projects` looks thin because the CSV has no study-level metadata. It still exists so samples can be grouped by study without parsing IDs, and so a later `start_date` or `indication` column has a place to go.

Treatment is subject-level **in this file**. If a person could receive miraclib and then something else, treatment would move to a visit table. That is not what the data shows today.

### How it would scale

SQLite with these indexes is enough for a take-home and for a few million `cell_counts` rows. For hundreds of projects:

- Keep the same four tables. Add `project_id` on `samples` only if you need to avoid a join on the hottest queries.
- Load CSV in chunks instead of one `read_csv`.
- Turn on WAL mode (`PRAGMA journal_mode=WAL`) for concurrent dashboard reads.
- If writes become a bottleneck, move to Postgres and keep the same schema; the SQL in `analysis.py` is ordinary parameterized SQL.

## Code structure

| File | What it does |
| --- | --- |
| `load_data.py` | Checks required columns, builds the schema, loads every CSV row |
| `analysis.py` | Frequencies, Mann-Whitney tests, boxplot, Part 4 SQL, B-cell average |
| `dashboard.py` | Reads the database plus `outputs/` and shows the results |
| `Makefile` | `setup`, `pipeline`, `dashboard` |

No extra packages, no class hierarchy. `analysis.py` is a list of functions in the order of the assignment.

## Statistical method

Question: in melanoma, miraclib, PBMC samples with a yes/no response, do the five relative frequencies differ between responders and non-responders?

Relative frequency is `count / (b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte) * 100`.

The test is two-sided Mann-Whitney U (`scipy.stats.mannwhitneyu`). Percentages are bounded and not assumed normal, so a rank test is the default. Five tests share one question, so p-values are Benjamini-Hochberg adjusted (`statsmodels.stats.multitest.multipletests`, `method="fdr_bh"`). A population is flagged significant when the adjusted p is below 0.05.

The assignment did not restrict to day 0, so **all matching samples** are used (three time points per subject). Those samples are not independent. The p-values should be read with that in mind; a stricter version would average each subject first.

## Actual results

Percentages per sample sum to 100% (max floating-point error `1.42e-14`).

### Responder vs non-responder (melanoma, miraclib, PBMC)

n = 993 responder samples, 975 non-responder samples.

| population | responder mean / median % | non-responder mean / median % | U statistic | p | BH p | significant |
| --- | --- | --- | --- | --- | --- | --- |
| b_cell | 9.80 / 9.43 | 10.00 / 9.79 | 459971 | 0.056 | 0.139 | no |
| cd8_t_cell | 24.88 / 24.73 | 24.94 / 24.60 | 478176 | 0.639 | 0.639 | no |
| cd4_t_cell | 30.54 / 30.22 | 29.90 / 29.66 | 515276 | 0.013 | 0.067 | no |
| nk_cell | 14.84 / 14.51 | 15.07 / 14.80 | 464546 | 0.121 | 0.202 | no |
| monocyte | 19.94 / 19.61 | 20.08 / 19.94 | 466510 | 0.163 | 0.204 | no |

CD4 T cells look slightly higher in responders on the raw p-value, but **nothing survives BH at 0.05**. The boxplot is `outputs/responder_boxplot.png`.

### Baseline subset (melanoma, PBMC, miraclib, day 0)

Queried from SQLite, not from the CSV in pandas.

Samples by project:

| project | n_samples |
| --- | --- |
| prj1 | 384 |
| prj3 | 272 |

prj2 does not appear because every prj2 sample in this file is whole blood (`WB`).

Unique subjects by response: no = 325, yes = 331.

Unique subjects by gender (`sex` in the CSV): F = 312, M = 344.

### B-cell average

Melanoma, male (`sex = M`), `response = yes`, `time_from_treatment_start = 0`, **all treatments and sample types** (not limited to miraclib or PBMC):

**10206.15** (n = 485 samples)

## Dashboard

After `make pipeline`:

```bash
make dashboard
```

The page shows dataset counts, the frequency table, the boxplot, the statistics table, the three baseline queries, and the B-cell average. It does not recompute the tests; it reads `outputs/` so what you see matches the CSV files.
