"""Load cell-count.csv into a small SQLite database.

Run from the repo root:
    python load_data.py
"""

from pathlib import Path

import pandas as pd
import sqlite3

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_count.db"

REQUIRED_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")

    df = pd.read_csv(path)
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"CSV is missing required columns: {missing}")

    if df["sample"].duplicated().any():
        raise ValueError("CSV has duplicate sample IDs; expected one row per sample.")

    for col in POPULATIONS:
        if (df[col] < 0).any():
            raise ValueError(f"Negative counts found in {col}")

    # Healthy donors have a blank response in the CSV. Keep that as NULL.
    df["response"] = df["response"].where(df["response"].isin(["yes", "no"]))
    return df


def create_schema(conn: sqlite3.Connection) -> None:
    # Drop first so reruns replace the whole database instead of appending rows.
    conn.executescript(
        """
        PRAGMA foreign_keys = OFF;
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;
        DROP TABLE IF EXISTS projects;
        PRAGMA foreign_keys = ON;

        CREATE TABLE projects (
            project_id TEXT PRIMARY KEY
        );

        CREATE TABLE subjects (
            subject_id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER NOT NULL,
            sex TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT,
            FOREIGN KEY (project_id) REFERENCES projects(project_id)
        );

        CREATE TABLE samples (
            sample_id TEXT PRIMARY KEY,
            subject_id TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL,
            FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
        );

        CREATE TABLE cell_counts (
            sample_id TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population),
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
        );

        CREATE INDEX idx_subjects_project ON subjects(project_id);
        CREATE INDEX idx_subjects_filters ON subjects(condition, treatment, response, sex);
        CREATE INDEX idx_samples_subject ON samples(subject_id);
        CREATE INDEX idx_samples_filters ON samples(sample_type, time_from_treatment_start);
        CREATE INDEX idx_counts_population ON cell_counts(population);
        """
    )


def insert_data(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    projects = pd.DataFrame({"project_id": sorted(df["project"].unique())})
    projects.to_sql("projects", conn, if_exists="append", index=False)

    subjects = (
        df[
            [
                "subject",
                "project",
                "condition",
                "age",
                "sex",
                "treatment",
                "response",
            ]
        ]
        .drop_duplicates("subject")
        .rename(columns={"subject": "subject_id", "project": "project_id"})
    )
    subjects.to_sql("subjects", conn, if_exists="append", index=False)

    samples = df[
        ["sample", "subject", "sample_type", "time_from_treatment_start"]
    ].rename(columns={"sample": "sample_id", "subject": "subject_id"})
    samples.to_sql("samples", conn, if_exists="append", index=False)

    count_frames = []
    for population in POPULATIONS:
        piece = df[["sample", population]].rename(
            columns={"sample": "sample_id", population: "count"}
        )
        piece["population"] = population
        count_frames.append(piece[["sample_id", "population", "count"]])
    cell_counts = pd.concat(count_frames, ignore_index=True)
    cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)


def main() -> None:
    df = load_csv(CSV_PATH)
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("PRAGMA foreign_keys = ON")
        create_schema(conn)
        insert_data(conn, df)
        conn.commit()

    print(f"Loaded {len(df)} samples into {DB_PATH.name}")
    print(f"  projects: {df['project'].nunique()}")
    print(f"  subjects: {df['subject'].nunique()}")
    print(f"  cell-count rows: {len(df) * len(POPULATIONS)}")


if __name__ == "__main__":
    main()
