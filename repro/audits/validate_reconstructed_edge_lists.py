#!/usr/bin/env python
"""Validate every canonical topology edge list in the review release."""
from __future__ import annotations

import argparse
import hashlib
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def canonical_edge_array(edge_index: np.ndarray) -> np.ndarray:
    edges = np.asarray(edge_index, dtype=np.int64)
    unique = np.unique(edges, axis=0)
    return unique[np.lexsort((unique[:, 1], unique[:, 0]))]


def edge_hash(edge_index: np.ndarray) -> str:
    return hashlib.sha256(canonical_edge_array(edge_index).astype("<i8", copy=False).tobytes()).hexdigest()


def degree_hash(edge_index: np.ndarray, n_nodes: int) -> str:
    edges = canonical_edge_array(edge_index)
    out = np.bincount(edges[:, 0], minlength=n_nodes)
    inc = np.bincount(edges[:, 1], minlength=n_nodes)
    return hashlib.sha256(np.concatenate([out, inc]).astype("<i8", copy=False).tobytes()).hexdigest()


def validate_root(root: Path) -> int:
    edge_root = root / "edge_lists" if (root / "edge_lists").is_dir() else root
    audit_root = edge_root.parent if edge_root.name == "edge_lists" else root
    audit_path = next(
        (
            candidate
            for candidate in (
                audit_root / "formal_edge_constraint_audit.csv",
                audit_root / "exact_edge_export_audit.csv",
                audit_root / "reconstructed_edge_audit.csv",
            )
            if candidate.exists()
        ),
        None,
    )
    if audit_path is None:
        raise FileNotFoundError(
            "no edge audit found; expected formal_edge_constraint_audit.csv, "
            "exact_edge_export_audit.csv or reconstructed_edge_audit.csv"
        )
    audit = pd.read_csv(audit_path)
    failures: list[str] = []
    edge_files = sorted((edge_root / "HMP2_to_PRISM").rglob("*.csv")) + sorted(
        (edge_root / "PRISM_to_HMP2").rglob("*.csv")
    )
    edge_files = [p for p in edge_files if p.name not in {"reconstructed_edge_audit.csv", "reconstructed_edge_summary.csv"}]
    if len(edge_files) != len(audit):
        failures.append(f"edge file count {len(edge_files)} != audit rows {len(audit)}")
    for _, row in audit.iterrows():
        direction = str(row.direction)
        split = int(row.split)
        condition = str(row.condition)
        if "edge_file" in audit.columns and str(row.edge_file) not in {"", "nan", "None"}:
            path = edge_root / str(row.edge_file)
        else:
            seed = "" if pd.isna(row.graph_seed) else f"__seed_{int(row.graph_seed)}"
            path = edge_root / direction / f"split_{split:02d}" / f"{condition}{seed}.csv"
        if not path.exists():
            failures.append(f"missing {path}")
            continue
        frame = pd.read_csv(path)
        if list(frame.columns) != ["source", "target"]:
            failures.append(f"bad columns in {path.name}")
            continue
        edges = frame[["source", "target"]].to_numpy(dtype=np.int64)
        n_nodes = int(row.n_nodes)
        unique = np.unique(edges, axis=0)
        self_loops = int(np.sum(edges[:, 0] == edges[:, 1]))
        out = np.bincount(edges[:, 0], minlength=n_nodes)
        if len(edges) != int(row.edge_count) or len(unique) != len(edges):
            failures.append(f"duplicate/count mismatch in {path.name}")
        if condition == "identity":
            if self_loops != n_nodes or len(edges) != n_nodes or not np.all(out == 1):
                failures.append(f"identity constraint failure in {path.name}")
        else:
            if self_loops != 0 or len(edges) != 8 * n_nodes or not np.all(out == 8):
                failures.append(f"non-identity constraint failure in {path.name}")
        if edge_hash(edges) != str(row.edge_hash) or degree_hash(edges, n_nodes) != str(row.degree_hash):
            failures.append(f"audit hash failure in {path.name}")
    print(f"edge_files={len(edge_files)} audit_rows={len(audit)}")
    if failures:
        print("FAIL:", "; ".join(failures))
        return 2
    print("PASS: all canonical edge lists satisfy local constraints and audit hashes")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--root", type=Path, help="extracted topology directory or its topology_constraints parent")
    group.add_argument("--archive", type=Path, help="formal_topology_edge_export_*.zip")
    args = parser.parse_args()
    if args.root is not None:
        return validate_root(args.root)
    with tempfile.TemporaryDirectory(prefix="ibd_topology_audit_") as tmp:
        with zipfile.ZipFile(args.archive) as archive:
            archive.extractall(tmp)
        extracted = Path(tmp) / "formal_topology_edge_export"
        if not extracted.is_dir():
            candidates = [p for p in Path(tmp).iterdir() if p.is_dir()]
            if len(candidates) != 1:
                raise FileNotFoundError("archive does not contain a unique topology export directory")
            extracted = candidates[0]
        return validate_root(extracted)


if __name__ == "__main__":
    raise SystemExit(main())
