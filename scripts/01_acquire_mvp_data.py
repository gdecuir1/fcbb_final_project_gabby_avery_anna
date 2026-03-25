"""
01_acquire_mvp_data.py

Phase 1 + Phase 2 (setup): Acquire the MVP cohort data.

Inputs
-----
- `config/cohort_mvp.yaml`
  - `cohort_study_id` (recommended: "TCGA-LUAD")
  - `paths.raw_dir` (where to store downloaded artifacts)
- One of:
  - cbioportal API configuration (preferred if available in your environment), OR
  - direct downloads of MAF files (fallback).

Outputs
--------
You should create (fill in exact filenames after you inspect available fields/endpoints):
- `data/raw/<cohort>_mutations.maf(.gz)?`
  - Raw mutation calls for primary solid tumors.
- `data/raw/<cohort>_sample_metadata.(tsv/csv)(.gz)?`
  - Sample identifiers and tumor type / primary status fields needed for filtering.

Functionality (to implement)
-------------------------------
- Use the cbioportal REST API (via `requests`) or download the needed MAFs.
- Store raw files into `data/raw/` without preprocessing.

Notes / things to watch
------------------------
- MAF schemas can vary by source; later scripts depend on which columns exist
  (e.g., variant class/type, coding region indicator, sample identifier key).
- Keep a small README/log of:
  - API query parameters or download URLs
  - date/time of download
  - which tumors were included/excluded.
"""

from __future__ import annotations

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Acquire MVP cohort mutation data.")
    parser.add_argument("--config", default="config/cohort_mvp.yaml", help="Path to cohort config YAML.")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    # TODO: Read YAML, query cbioportal (or download MAF), write to data/raw/.
    raise NotImplementedError("Fill in cbioportal acquisition / MAF download logic.")


if __name__ == "__main__":
    main()

