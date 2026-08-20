"""Streamlit dashboard for the cell-count analysis.

Run from the repo root:
    streamlit run dashboard.py
"""

from pathlib import Path

import pandas as pd
import sqlite3
import streamlit as st

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_count.db"
OUTPUT_DIR = ROOT / "outputs"


def require_outputs() -> None:
    needed = [
        OUTPUT_DIR / "cell_frequencies.csv",
        OUTPUT_DIR / "statistical_results.csv",
        OUTPUT_DIR / "responder_boxplot.png",
        OUTPUT_DIR / "baseline_samples_by_project.csv",
        OUTPUT_DIR / "baseline_subjects_by_response.csv",
        OUTPUT_DIR / "baseline_subjects_by_gender.csv",
        OUTPUT_DIR / "b_cell_average.csv",
    ]
    missing = [p.name for p in needed if not p.exists()]
    if missing:
        st.error(
            "Output files are missing: "
            + ", ".join(missing)
            + ". Run `make pipeline` first."
        )
        st.stop()


def dataset_overview() -> dict:
    conn = sqlite3.connect(DB_PATH)
    overview = {
        "projects": conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0],
        "subjects": conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0],
        "samples": conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0],
        "cell_count_rows": conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0],
    }
    by_condition = pd.read_sql_query(
        "select condition, count(*) AS n_subjects from subjects group by condition order by condition",
        conn,
    )
    by_treatment = pd.read_sql_query(
        "select treatment, count(*) AS n_subjects from subjects group by treatment order by treatment",
        conn,
    )
    by_sample_type = pd.read_sql_query(
        "select sample_type, count(*) AS n_samples from samples group by sample_type order by sample_type",
        conn,
    )
    conn.close()
    return overview, by_condition, by_treatment, by_sample_type


st.set_page_config(page_title="Cell-count analysis", layout="wide")
st.title("Immune cell-count analysis")
st.caption("Melanoma / carcinoma PBMC and whole-blood counts loaded from cell-count.csv")

require_outputs()
overview, by_condition, by_treatment, by_sample_type = dataset_overview()

st.header("Dataset overview")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Projects", overview["projects"])
c2.metric("Subjects", overview["subjects"])
c3.metric("Samples", overview["samples"])
c4.metric("Cell-count rows", overview["cell_count_rows"])

left, mid, right = st.columns(3)
with left:
    st.subheader("Subjects by condition")
    st.dataframe(by_condition, hide_index=True)
with mid:
    st.subheader("Subjects by treatment")
    st.dataframe(by_treatment, hide_index=True)
with right:
    st.subheader("Samples by type")
    st.dataframe(by_sample_type, hide_index=True)

st.header("Cell frequency table")
st.write(
    "Long-format percentages for every sample: "
    "`percentage = count / (b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte) * 100`."
)
freq = pd.read_csv(OUTPUT_DIR / "cell_frequencies.csv")
st.dataframe(freq.head(100), hide_index=True)
st.caption(f"Showing the first 100 of {len(freq):,} rows. Full table: outputs/cell_frequencies.csv")

st.header("Responder vs non-responder boxplot")
st.write("Relative frequencies for melanoma + miraclib + PBMC samples, split by response.")
st.image(str(OUTPUT_DIR / "responder_boxplot.png"), use_container_width=True)

st.header("Statistical results")
st.write(
    "Two-sided Mann-Whitney U on relative frequencies, with Benjamini-Hochberg "
    "adjusted p-values. Significant if adjusted p < 0.05."
)
stats = pd.read_csv(OUTPUT_DIR / "statistical_results.csv")
st.dataframe(stats, hide_index=True)

st.header("Baseline subset analysis")
st.write("SQLite query: melanoma, PBMC, miraclib, time_from_treatment_start = 0.")
b1, b2, b3 = st.columns(3)
with b1:
    st.subheader("Samples by project")
    st.dataframe(pd.read_csv(OUTPUT_DIR / "baseline_samples_by_project.csv"), hide_index=True)
with b2:
    st.subheader("Subjects by response")
    st.dataframe(pd.read_csv(OUTPUT_DIR / "baseline_subjects_by_response.csv"), hide_index=True)
with b3:
    st.subheader("Subjects by gender")
    st.dataframe(pd.read_csv(OUTPUT_DIR / "baseline_subjects_by_gender.csv"), hide_index=True)
st.caption("Gender is the CSV `sex` column (M/F). Project prj2 is whole blood only, so it does not appear here.")

st.header("Final B-cell average")
st.write(
    "Melanoma males, response = yes, time_from_treatment_start = 0, "
    "across all treatments and sample types. Not filtered to miraclib or PBMC."
)
bcell = pd.read_csv(OUTPUT_DIR / "b_cell_average.csv")
st.metric("Average b_cell count", f"{bcell.loc[0, 'mean_b_cell']:.2f}")
st.caption(f"n = {int(bcell.loc[0, 'n_samples'])} samples")
