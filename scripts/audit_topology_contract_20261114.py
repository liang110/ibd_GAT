#!/usr/bin/env python3
"""Audit exported directed graph edge lists against the topology contract.

The script is deliberately independent of model predictions.  It checks the
edge objects that are sufficient to verify non-self simple graphs, fixed
out-degree and, when rewired edge lists are supplied, preservation of the
directed in/out-degree sequence relative to the paired real graph.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _read_edges(path: Path) -> list[tuple[int, int]]:
    frame = pd.read_csv(path)
    required = {"source", "target"}
    if not required.issubset(frame.columns):
        raise ValueError(f"{path}: expected columns source,target")
    return [(int(s), int(t)) for s, t in frame[["source", "target"]].itertuples(index=False)]


def _degree(edges: list[tuple[int, int]], n_nodes: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    out = [0] * n_nodes
    inc = [0] * n_nodes
    for source, target in edges:
        out[source] += 1
        inc[target] += 1
    return tuple(out), tuple(inc)


def _audit(path: Path, expected_k: int | None) -> dict[str, object]:
    edges = _read_edges(path)
    if not edges:
        raise ValueError(f"{path}: empty edge list")
    nodes = sorted({node for edge in edges for node in edge})
    n_nodes = max(nodes) + 1
    unique = len(set(edges)) == len(edges)
    self_loops = sum(source == target for source, target in edges)
    out, inc = _degree(edges, n_nodes)
    if expected_k is not None and any(value != expected_k for value in out):
        raise ValueError(f"{path}: out-degree is not fixed at {expected_k}: {sorted(set(out))}")
    if self_loops:
        raise ValueError(f"{path}: contains {self_loops} self-loop(s)")
    if not unique:
        raise ValueError(f"{path}: contains duplicate directed edges")
    return {
        "file": path.name,
        "n_nodes": n_nodes,
        "edge_count": len(edges),
        "unique_directed_edges": unique,
        "self_loop_count": self_loops,
        "out_degree_min": min(out),
        "out_degree_max": max(out),
        "in_degree_min": min(inc),
        "in_degree_max": max(inc),
        "out_degree": out,
        "in_degree": inc,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-k", type=int, default=8)
    parser.add_argument(
        "--require-rewired",
        action="store_true",
        help="fail if no files whose name contains 'rewired' are supplied",
    )
    args = parser.parse_args()

    files = sorted(args.edges_root.glob("*_edges.csv"))
    if not files:
        raise SystemExit(f"no *_edges.csv files found under {args.edges_root}")
    records = [_audit(path, args.expected_k) for path in files]
    rewired = [record for record in records if "rewired" in str(record["file"]).lower()]
    if args.require_rewired and not rewired:
        raise SystemExit("rewired edge lists are required but none were supplied")

    paired_checks: list[dict[str, object]] = []
    by_key = {record["file"]: record for record in records}
    for real_name, real in by_key.items():
        if "real" not in real_name.lower() and "spearman" not in real_name.lower():
            continue
        stem = real_name.replace("_spearman_edges.csv", "").replace("_real_edges.csv", "")
        for rewired_name, null in by_key.items():
            if "rewired" not in rewired_name.lower() or not rewired_name.startswith(stem):
                continue
            paired_checks.append(
                {
                    "real_file": real_name,
                    "rewired_file": rewired_name,
                    "same_edge_count": real["edge_count"] == null["edge_count"],
                    "same_out_degree": real["out_degree"] == null["out_degree"],
                    "same_in_degree": real["in_degree"] == null["in_degree"],
                    "different_edge_set": real.get("file") != null.get("file"),
                }
            )

    result = {
        "n_files": len(records),
        "n_rewired_files": len(rewired),
        "expected_k": args.expected_k,
        "records": records,
        "paired_rewired_checks": paired_checks,
        "rewired_edge_lists_present": bool(rewired),
        "status": "PASS" if (rewired or not args.require_rewired) else "INCOMPLETE",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("n_files", "n_rewired_files", "status")}, indent=2))


if __name__ == "__main__":
    main()
