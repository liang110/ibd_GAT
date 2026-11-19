"""Compare the synthetic example output with the expected paired scores."""

import argparse
from pathlib import Path

import pandas as pd


def load_scores(path):
    keys = ["direction", "target_stratum", "scoring_unit", "candidate", "reference", "metric"]
    return pd.read_csv(path).set_index(keys, verify_integrity=True).sort_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="CSV produced by score_fixed_predictions.py")
    args = parser.parse_args()
    expected = load_scores(Path(__file__).with_name("expected_scores.csv"))
    actual = load_scores(args.output)
    pd.testing.assert_frame_equal(actual, expected, check_exact=False, rtol=0, atol=1e-12)
    print(f"Verified all {len(expected)} comparisons.")


if __name__ == "__main__":
    main()
