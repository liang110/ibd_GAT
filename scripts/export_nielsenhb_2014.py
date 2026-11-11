#!/usr/bin/env python
"""Export the pinned public NielsenHB_2014 resources without changing them."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ibd_microbiome_framework.data.nielsenhb_2014 import export_public_resources


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relative-abundance-rda", type=Path, required=True)
    parser.add_argument("--sample-metadata-rda", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = export_public_resources(
        args.relative_abundance_rda.resolve(),
        args.sample_metadata_rda.resolve(),
        args.output_dir.resolve(),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
