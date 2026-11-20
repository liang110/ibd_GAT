"""Verify the fitted graph inputs, their degree controls and taxon-aligned overlap."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


DIRECTIONS = ("HMP2_to_PRISM", "PRISM_to_HMP2")
FIT_KEYS = ["direction", "split"]
GRAPH_KEYS = FIT_KEYS + ["condition", "graph_seed"]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def int_digest(values):
    return hashlib.sha256(np.asarray(values, dtype="<i8").tobytes()).hexdigest()


def verify_graph(row, edge_array, names):
    key = (row.direction, row.split, row.condition, row.graph_seed)
    count = len(names)
    require(count == row.n_nodes and len(edge_array) == row.n_edges, f"Graph size: {key}")
    require(np.all((edge_array >= 0) & (edge_array < count)), f"Node outside fit: {key}")
    ordered = np.unique(edge_array, axis=0)
    require(len(ordered) == len(edge_array), f"Duplicate edge: {key}")
    src, tgt = ordered.T
    out_degree = np.bincount(src, minlength=count)
    in_degree = np.bincount(tgt, minlength=count)
    if row.condition == "identity":
        require(np.array_equal(src, tgt) and np.all(out_degree == 1), f"Self-loop control: {key}")
    else:
        require(np.all(src != tgt) and np.all(out_degree == 8), f"Non-self graph: {key}")
    feature_hash = hashlib.sha256("\n".join(names).encode("utf8")).hexdigest()
    require(feature_hash == row.checkpoint_feature_hash == row.selected_feature_hash, f"Taxon order: {key}")
    require(int_digest(ordered) == row.edge_hash == row.edge_hash_fit, f"Fit edge identity: {key}")
    require(int_digest(np.r_[out_degree, in_degree]) == row.degree_hash == row.degree_hash_fit, f"Degree sequence: {key}")
    return frozenset(zip(src, tgt))


def verify_overlap(real_graphs, taxon_order, summary):
    for direction in DIRECTIONS:
        local = [real_graphs[(direction, split)] for split in range(1, 11)]
        mapped = []
        for split, edges in enumerate(local, 1):
            names = taxon_order[(direction, split)]
            mapped.append({(names[source], names[target]) for source, target in edges})
        row = summary.loc[direction]
        for column, graphs in [("local_mean", local), ("taxon_mean", mapped)]:
            values = [len(left & right) / len(left | right) for left, right in itertools.combinations(graphs, 2)]
            require(len(values) == row.n_pairs == 45 and abs(np.mean(values) - row[column]) < 1e-12,
                    f"Between-fit overlap: {direction}, {column}")
            if column == "taxon_mean":
                require(abs(min(values) - row.taxon_min) < 1e-12 and abs(max(values) - row.taxon_max) < 1e-12,
                        f"Taxon overlap range: {direction}")


def validate_topology(root):
    path = root / "results/topology"
    features = pd.read_csv(path / "features.csv")
    expected = {(direction, split) for direction in DIRECTIONS for split in range(1, 11)}
    names = {}
    for key, group in features.groupby(FIT_KEYS):
        ordered = group.sort_values("node")
        require(ordered.node.tolist() == list(range(len(ordered))) and ordered.taxon.is_unique,
                f"Node ordering: {key}")
        names[key] = ordered.taxon.tolist()
    require(set(names) == expected, "Expected two directions and ten feature sets each")
    metrics_file = path / "edge_ablation_metrics.csv"
    training = json.loads((path / "training.json").read_text(encoding="utf8"))
    normalized = metrics_file.read_bytes().replace(b"\r\n", b"\n")
    require(hashlib.sha256(normalized).hexdigest() == training["fit_metrics_sha256_lf"], "Fitted metrics checksum")
    fits = pd.read_csv(metrics_file)
    index = pd.read_csv(path / "graphs.csv")
    graphs = index.merge(fits, on=GRAPH_KEYS, how="outer", validate="one_to_one", suffixes=("", "_fit"), indicator=True)
    require(len(graphs) == 140 and graphs._merge.eq("both").all() and graphs.graph_id.is_unique,
            "Unmatched graph or fitted model")
    require(set(graphs[FIT_KEYS].itertuples(index=False,name=None)) == expected, "Missing source fit")
    edge_table = pd.read_csv(path / "edges.csv", dtype={"graph_id":"int64", "source":"int64", "target":"int64"})
    arrays = {gid:g[["source","target"]].to_numpy() for gid,g in edge_table.groupby("graph_id")}
    require(set(graphs.graph_id) == set(arrays), "Edge rows and graph IDs differ")
    real_graphs = {}
    for key, group in graphs.groupby(FIT_KEYS):
        require(group.condition.value_counts().to_dict() == {"rewired":5, "identity":1, "real_spearman":1},
                f"Graph conditions: {key}")
        real = group.loc[group.condition == "real_spearman"].iloc[0]
        rewired = group.loc[group.condition == "rewired"]
        require(rewired.graph_seed.notna().all() and rewired.graph_seed.nunique() == 5 and
                rewired.degree_hash.eq(real.degree_hash).all() and rewired.edge_hash.ne(real.edge_hash).all(),
                f"Degree-preserving rewiring: {key}")
        for row in group.itertuples(index=False):
            directed = verify_graph(row, arrays[row.graph_id], names[key])
            if row.condition == "real_spearman":
                real_graphs[key] = directed
    summary = pd.read_csv(root / "results/comparisons/topology_taxon_aligned_summary.csv").set_index("direction",verify_integrity=True)
    verify_overlap(real_graphs,names,summary)
    print(f"Verified {len(graphs)} graph inputs ({len(edge_table):,} directed edges), 20 taxon orders and 45 fit pairs per direction.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    validate_topology(parser.parse_args().root)
