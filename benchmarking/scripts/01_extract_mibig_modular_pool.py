#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

TARGET_CLASSES = {"NRP", "Polyketide"}
PRIOR_MIBIG = {
    "BGC0002072",
    "BGC0002068",
    "BGC0000320",
    "BGC0002124",
    "BGC0001056",
    "BGC0000921",
    "BGC0002082",
}


def module_count(cluster: dict) -> tuple[int, int]:
    nrp_count = 0
    nrp = cluster.get("nrp")
    if isinstance(nrp, dict):
        genes = nrp.get("nrps_genes", [])
        if isinstance(genes, list):
            for gene in genes:
                if isinstance(gene, dict):
                    modules = gene.get("modules", [])
                    if isinstance(modules, list):
                        nrp_count += len(modules)

    pks_count = 0
    polyketide = cluster.get("polyketide")
    if isinstance(polyketide, dict):
        synthases = polyketide.get("synthases", [])
        if isinstance(synthases, list):
            for synthase in synthases:
                if isinstance(synthase, dict):
                    modules = synthase.get("modules", [])
                    if isinstance(modules, list):
                        pks_count += len(modules)

    return nrp_count, pks_count


def compound_names(cluster: dict) -> str:
    values = []
    compounds = cluster.get("compounds", [])
    if isinstance(compounds, list):
        for item in compounds:
            if isinstance(item, dict):
                name = str(item.get("compound", "")).strip()
                if name and name not in values:
                    values.append(name)
    return "; ".join(values)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract active modular NRPS/PKS MIBiG entries with genomic loci "
            "for prospective BGCPhaser validation discovery."
        )
    )
    parser.add_argument("--mibig-data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--min-modules",
        type=int,
        default=3,
        help="Minimum annotated NRP+PKS module count; default 3.",
    )
    args = parser.parse_args()

    if args.min_modules < 1:
        raise SystemExit("--min-modules must be >= 1")
    if not args.mibig_data.is_dir():
        raise SystemExit(f"MIBiG data directory absent: {args.mibig_data}")

    rows = []
    for path in sorted(args.mibig_data.glob("BGC*.json")):
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        cluster = data.get("cluster")
        if not isinstance(cluster, dict):
            continue
        if cluster.get("status") != "active":
            continue

        accession = str(cluster.get("mibig_accession", "")).strip()
        if not accession or accession in PRIOR_MIBIG:
            continue

        classes = cluster.get("biosyn_class", [])
        if not isinstance(classes, list):
            continue
        class_set = {str(value) for value in classes}
        if not (class_set & TARGET_CLASSES):
            continue

        loci = cluster.get("loci")
        if not isinstance(loci, dict):
            continue
        reference = str(loci.get("accession", "")).strip()
        start = loci.get("start_coord")
        end = loci.get("end_coord")
        if not reference or not isinstance(start, int) or not isinstance(end, int):
            continue
        if end <= start:
            continue

        tax_id = str(cluster.get("ncbi_tax_id", "")).strip()
        organism = str(cluster.get("organism_name", "")).strip()
        if not tax_id or not organism:
            continue

        nrp_modules, pks_modules = module_count(cluster)
        total_modules = nrp_modules + pks_modules
        if total_modules < args.min_modules:
            continue

        completeness = str(loci.get("completeness", "")).strip() or "NA"
        rows.append(
            {
                "mibig_accession": accession,
                "organism_name": organism,
                "ncbi_tax_id": tax_id,
                "biosyn_class": ";".join(sorted(class_set)),
                "reference_accession": reference,
                "start_coord": start,
                "end_coord": end,
                "locus_span_bp": end - start,
                "locus_completeness": completeness,
                "nrp_module_count": nrp_modules,
                "pks_module_count": pks_modules,
                "total_module_count": total_modules,
                "compounds": compound_names(cluster),
            }
        )

    rows.sort(
        key=lambda row: (
            -int(row["total_module_count"]),
            -int(row["locus_span_bp"]),
            str(row["mibig_accession"]),
        )
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "mibig_accession",
        "organism_name",
        "ncbi_tax_id",
        "biosyn_class",
        "reference_accession",
        "start_coord",
        "end_coord",
        "locus_span_bp",
        "locus_completeness",
        "nrp_module_count",
        "pks_module_count",
        "total_module_count",
        "compounds",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} modular MIBiG candidates to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
