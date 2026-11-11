#!/usr/bin/env python
"""Create immutable data-contract artifacts from the frozen legacy inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ibd_microbiome_framework.data.contracts import MatrixSpec, write_contracts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.legacy_root.resolve()
    specs = [
        MatrixSpec("HMP2", "microbiome", root / "manual_abundance.csv", "Feature", root / "clean_metadata.csv"),
        MatrixSpec("HMP2", "metabolome", root / "metabolites.csv", "Feature", root / "clean_metadata.csv"),
        MatrixSpec("PRISM", "microbiome", root / "prism_abundance_clean.csv", "Species", root / "prism_metadata_clean.csv"),
        MatrixSpec("PRISM", "metabolome", root / "prism_metabolites_clean.csv", "Metabolite", root / "prism_metadata_clean.csv"),
    ]
    summary = write_contracts(specs, args.output_dir.resolve(), root)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
