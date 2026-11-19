"""Rebuild the two microbiome panels from their distributed profiles and metadata."""
import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

from ibd_microbiome_framework.data.contracts import MatrixSpec, write_contracts
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list


def hmp_inputs(directory):
    meta = pd.read_csv(directory / "hmp2_metadata_2018-08-20.csv", low_memory=False)
    meta = meta.loc[(meta["data_type"] == "metagenomics") & meta["External ID"].notna()].copy()
    meta = meta.rename(columns={"External ID": "sample_id"})
    meta["group"] = meta["diagnosis"].map({"CD": "CD", "UC": "UC", "nonIBD": "Control"})
    if meta["group"].isna().any() or meta["sample_id"].duplicated().any():
        raise ValueError("HMP2 diagnoses or External IDs do not match the published MGX release")
    taxa = pd.read_csv(directory / "taxonomic_profiles_3.tsv.gz", sep="\t", index_col=0, low_memory=False)
    taxa.columns = taxa.columns.str.removesuffix("_profile")
    species = taxa.index.str.contains("s__", regex=False) & ~taxa.index.str.contains("t__", regex=False)
    taxa = taxa.loc[species].copy()
    taxa.index = taxa.index.str.split("|").str[-1]
    columns = taxa.columns.intersection(meta["sample_id"], sort=True)
    if taxa.columns.duplicated().any() or len(columns) != 1638 or len(taxa) != 578:
        raise ValueError("Expected 578 HMP2 reported species rows and 1,638 matched MGX records")
    return taxa[columns], meta.loc[meta.sample_id.isin(columns)]


def franzosa_inputs(directory):
    path = directory / "41564_2018_306_MOESM6_ESM.xlsx"
    with path.open("rb") as handle:
        workbook = load_workbook(handle, read_only=True, data_only=True)
        try:
            sheet = workbook["dataset_s4"]
            rows = iter(sheet.values)
            next(rows)  # Dataset title; the second row contains the study sample IDs.
            columns = next(rows)
            profiles = pd.DataFrame(rows, columns=columns).set_index(columns[0])
        finally:
            workbook.close()
    profiles.index.name = "Species"
    profiles = profiles.replace("#N/A", float("nan"))
    diagnosis = profiles.loc["Diagnosis"]
    if profiles.shape != (209, 220) or not diagnosis.isin(["CD", "UC", "Control"]).all():
        raise ValueError("Unexpected Franzosa Dataset 4 dimensions or diagnoses")
    metadata = pd.DataFrame({"sample_id": profiles.columns, "group": diagnosis.to_numpy()})
    if metadata.sample_id.duplicated().any():
        raise ValueError("Duplicated Franzosa participant ID")
    return profiles, metadata


def prepare(raw, work, scored_ids):
    hmp, hmeta = hmp_inputs(raw)
    prism, pmeta = franzosa_inputs(raw)
    hmp.to_csv(work / "manual_abundance.csv")
    columns = ["sample_id", "Participant ID", "week_num", "visit_num", "group", "diagnosis", "site_name"]
    hmeta[columns].to_csv(work / "clean_metadata.csv", index=False)
    prism.to_csv(work / "prism_abundance_clean.csv")
    pmeta.to_csv(work / "prism_metadata_clean.csv", index=False)
    contract = work / "contract"
    specs = [
        MatrixSpec("HMP2", "microbiome", work / "manual_abundance.csv", "Feature", work / "clean_metadata.csv"),
        MatrixSpec("PRISM", "microbiome", work / "prism_abundance_clean.csv", "Species", work / "prism_metadata_clean.csv"),
    ]
    write_contracts(specs, contract, work)
    shared = load_shared_feature_list(contract)
    if hashlib.sha256("\n".join(shared).encode()).hexdigest() != "f1a8cda85f11fb9094743b49924b6ba93515c1f0cc09e947f76d60711b0f23b4":
        raise ValueError("Shared taxon identities or ordering differ from the study")
    panels = {}
    expected_profiles = {
        "HMP2": "ef67df595bfb12509a4f63f830027e6691f2b0dec68c3c1c1e08dc8170083056",
        "PRISM": "15fd8a6e23662180ac18734e62ce16dee82913d17265c6a8dc449e2ec2be247f",
    }
    for spec in specs:
        manifest = pd.read_csv(contract / f"cohort_manifest_{spec.cohort.lower()}.csv")
        panels[spec.cohort] = load_canonical_microbiome(spec.path, manifest, shared)
        serial = panels[spec.cohort].sort_index().to_csv(float_format="%.12f", lineterminator="\n")
        if hashlib.sha256(serial.encode("utf-8")).hexdigest() != expected_profiles[spec.cohort]:
            raise ValueError(f"{spec.cohort}: sample IDs or projected abundances differ from the frozen panel")
    retained = set(panels["HMP2"].index)
    retained_meta = hmeta[hmeta["sample_id"].isin(retained)]
    observed = [len(retained), retained_meta["Participant ID"].nunique(),
                retained_meta[["Participant ID", "visit_num"]].drop_duplicates().shape[0], len(panels["PRISM"])]
    if observed != [1626, 130, 1584, 220]:
        raise ValueError(f"Target panel check failed: {observed}")
    depth = hmeta.set_index("sample_id")["reads_filtered"]
    if depth.isna().any():
        raise ValueError("Missing HMP2 read-depth metadata")
    high = set(depth[depth >= 1000000].index)
    if len(high) != 1595 or not high <= retained or len(retained - high) != 31:
        raise ValueError("HMP2 read-depth and shared-feature inclusion sets differ from the study")
    if scored_ids:
        scored = pd.read_csv(scored_ids, sep="\t", header=None, usecols=[0])[0]
        if len(scored) != 1595 or set(scored) != high:
            raise ValueError("Public dysbiosis sample IDs differ from the 1,595 matched high-depth records")
    # Preserve relative paths after moving this staging directory into place.
    inventory = pd.read_csv(contract / "source_inventory.csv")
    for key in ["path", "metadata_path"]:
        inventory[key] = inventory[key].map(lambda x: str(Path(x).relative_to(work).as_posix()))
    inventory.to_csv(contract / "source_inventory.csv", index=False)
    summary_path = contract / "contract_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["source_root"] = ".. (relative to the contract directory)"
    summary["input_modalities"] = ["microbiome"]
    summary["shared_metabolite_exact_name_candidates"] = None
    summary["hmp2_included_records"] = len(retained)
    summary["hmp2_read_depth_at_least_one_million"] = len(high)
    summary["hmp2_included_below_one_million"] = len(retained - high)
    summary["hmp2_published_dysbiosis_ids_verified"] = scored_ids is not None
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory with the three published files.")
    parser.add_argument("--output-dir", required=True, type=Path, help="A new directory for inputs and contract/.")
    parser.add_argument("--published-dysbiosis-ids", type=Path, help="Optional headerless public sample-score TSV.")
    args = parser.parse_args()
    destination = args.output_dir.resolve()
    if destination.exists():
        parser.error("Choose a new output directory.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        work = Path(temporary) / "prepared"
        work.mkdir()
        prepare(args.input_dir, work, args.published_dysbiosis_ids)
        work.rename(destination)
    print("Prepared 126 shared candidate taxa, 1,626 HMP2 records and 220 Franzosa records.")


if __name__ == "__main__":
    main()
