"""Run the analysis steps against cell_count.db and write files under outputs/."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import sqlite3
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_count.db"
OUTPUT_DIR = ROOT / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def connect() -> sqlite3.Connection:
    if not DB_PATH.exists():
        raise FileNotFoundError("cell_count.db not found. Run python load_data.py first.")
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def sample_counts(conn: sqlite3.Connection) -> pd.DataFrame:
    """One row per sample, with the five population counts as columns."""
    long = pd.read_sql_query(
        "SELECT sample_id AS sample, population, count FROM cell_counts",
        conn,
    )
    wide = long.pivot(index="sample", columns="population", values="count")
    missing = [col for col in POPULATIONS if col not in wide.columns]
    if missing:
        raise ValueError(f"Database is missing cell populations: {missing}")
    return wide[POPULATIONS].reset_index()


def cell_frequencies(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.copy()
    counts["total_count"] = counts[POPULATIONS].sum(axis=1)
    if (counts["total_count"] <= 0).any():
        raise ValueError("Found a sample with total_count <= 0")

    long = counts.melt(
        id_vars=["sample", "total_count"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["percentage"] = long["count"] / long["total_count"] * 100

    checksum = long.groupby("sample")["percentage"].sum()
    max_error = (checksum - 100).abs().max()
    if max_error > 1e-6:
        raise ValueError(f"Percentages do not sum to 100% (max abs error={max_error})")
    print(f"Frequency check: percentages sum to 100% (max abs error={max_error:.2e})")
    return long.sort_values(["sample", "population"]).reset_index(drop=True)


def stats_subset(conn: sqlite3.Connection, frequencies: pd.DataFrame) -> pd.DataFrame:
    samples = pd.read_sql_query(
        """
        SELECT
            s.sample_id AS sample,
            u.response
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = ?
          AND u.treatment = ?
          AND s.sample_type = ?
          AND u.response IN (?, ?)
        """,
        conn,
        params=("melanoma", "miraclib", "PBMC", "yes", "no"),
    )
    merged = frequencies.merge(samples, on="sample", how="inner")
    return merged


def mann_whitney_table(subset: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population in POPULATIONS:
        pop = subset[subset["population"] == population]
        yes = pop.loc[pop["response"] == "yes", "percentage"]
        no = pop.loc[pop["response"] == "no", "percentage"]
        result = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append(
            {
                "population": population,
                "responder_count": int(yes.count()),
                "non_responder_count": int(no.count()),
                "responder_mean": yes.mean(),
                "responder_median": yes.median(),
                "non_responder_mean": no.mean(),
                "non_responder_median": no.median(),
                "statistic": result.statistic,
                "p_value": result.pvalue,
            }
        )

    table = pd.DataFrame(rows)
    _, adj, _, _ = multipletests(table["p_value"], method="fdr_bh")
    table["p_value_adj"] = adj
    table["significant"] = table["p_value_adj"] < 0.05
    return table


def save_boxplot(subset: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(14, 4.5), sharey=True)
    for ax, population in zip(axes, POPULATIONS):
        pop = subset[subset["population"] == population]
        no = pop.loc[pop["response"] == "no", "percentage"]
        yes = pop.loc[pop["response"] == "yes", "percentage"]
        ax.boxplot(
            [no, yes],
            tick_labels=["no", "yes"],
            widths=0.55,
            patch_artist=True,
            boxprops={"facecolor": "#d9e8f5"},
            medianprops={"color": "#1f4e79"},
        )
        ax.set_title(population)
        ax.set_xlabel("response")
        ax.grid(axis="y", linestyle=":", alpha=0.4)

    axes[0].set_ylabel("relative frequency (%)")
    fig.suptitle(
        "Responder vs non-responder cell frequencies\n"
        "melanoma, miraclib, PBMC",
        y=1.05,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def part4_queries(conn: sqlite3.Connection) -> dict[str, pd.DataFrame]:
    """Baseline subset: melanoma + PBMC + miraclib + day 0. Query SQLite directly."""
    filters = {
        "condition": "melanoma",
        "sample_type": "PBMC",
        "treatment": "miraclib",
        "time_from_treatment_start": 0,
    }

    samples_by_project = pd.read_sql_query(
        """
        SELECT u.project_id AS project, COUNT(s.sample_id) AS n_samples
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = :condition
          AND s.sample_type = :sample_type
          AND u.treatment = :treatment
          AND s.time_from_treatment_start = :time_from_treatment_start
        GROUP BY u.project_id
        ORDER BY u.project_id
        """,
        conn,
        params=filters,
    )

    subjects_by_response = pd.read_sql_query(
        """
        SELECT u.response, COUNT(DISTINCT u.subject_id) AS n_subjects
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = :condition
          AND s.sample_type = :sample_type
          AND u.treatment = :treatment
          AND s.time_from_treatment_start = :time_from_treatment_start
        GROUP BY u.response
        ORDER BY u.response
        """,
        conn,
        params=filters,
    )

    subjects_by_gender = pd.read_sql_query(
        """
        SELECT u.sex AS gender, COUNT(DISTINCT u.subject_id) AS n_subjects
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = :condition
          AND s.sample_type = :sample_type
          AND u.treatment = :treatment
          AND s.time_from_treatment_start = :time_from_treatment_start
        GROUP BY u.sex
        ORDER BY u.sex
        """,
        conn,
        params=filters,
    )
    return {
        "baseline_samples_by_project": samples_by_project,
        "baseline_subjects_by_response": subjects_by_response,
        "baseline_subjects_by_gender": subjects_by_gender,
    }


def b_cell_average(conn: sqlite3.Connection) -> pd.DataFrame:
    """Melanoma males, response=yes, day 0, all treatments and sample types."""
    result = pd.read_sql_query(
        """
        SELECT
            COUNT(*) AS n_samples,
            ROUND(AVG(c.count), 2) AS mean_b_cell
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        JOIN cell_counts AS c ON c.sample_id = s.sample_id
        WHERE u.condition = ?
          AND u.sex = ?
          AND u.response = ?
          AND s.time_from_treatment_start = ?
          AND c.population = ?
        """,
        conn,
        params=("melanoma", "M", "yes", 0, "b_cell"),
    )
    return result


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    conn = connect()

    counts = sample_counts(conn)
    frequencies = cell_frequencies(counts)
    frequencies.to_csv(OUTPUT_DIR / "cell_frequencies.csv", index=False)
    print(f"Wrote outputs/cell_frequencies.csv ({len(frequencies)} rows)")

    subset = stats_subset(conn, frequencies)
    stats = mann_whitney_table(subset)
    stats.to_csv(OUTPUT_DIR / "statistical_results.csv", index=False)
    print(f"Wrote outputs/statistical_results.csv ({len(stats)} populations)")

    save_boxplot(subset, OUTPUT_DIR / "responder_boxplot.png")
    print("Wrote outputs/responder_boxplot.png")

    for name, table in part4_queries(conn).items():
        path = OUTPUT_DIR / f"{name}.csv"
        table.to_csv(path, index=False)
        print(f"Wrote {path.relative_to(ROOT)}")
        print(table.to_string(index=False))

    bcell = b_cell_average(conn)
    bcell.to_csv(OUTPUT_DIR / "b_cell_average.csv", index=False)
    print("Wrote outputs/b_cell_average.csv")
    print(
        "Melanoma males, response=yes, day 0, all treatments/sample types: "
        f"mean b_cell = {bcell.loc[0, 'mean_b_cell']:.2f} "
        f"(n={int(bcell.loc[0, 'n_samples'])})"
    )

    conn.close()


if __name__ == "__main__":
    main()
