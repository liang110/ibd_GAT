#!/usr/bin/env python
"""Prepare NielsenHB_2014 only as a locked, source-feature external target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ibd_microbiome_framework.data.nielsenhb_2014 import prepare_locked_external_target


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relative-abundance-csv", type=Path, required=True)
    parser.add_argument("--sample-metadata-csv", type=Path, required=True)
    parser.add_argument("--source-shared-features-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = prepare_locked_external_target(
        args.relative_abundance_csv.resolve(),
        args.sample_metadata_csv.resolve(),
        args.source_shared_features_csv.resolve(),
        args.output_dir.resolve(),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
