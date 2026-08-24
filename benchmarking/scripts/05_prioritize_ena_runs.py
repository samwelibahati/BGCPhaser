#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


def norm_name(value: str) -> str:
    text = value.strip().lower()
    text = re.sub(r"[\[\](),]", " ", text)
    text = re.sub(r"\b(subsp\.?|ssp\.?|strain|str\.?)\b", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def parse_int(value: str) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def has_two_fastqs(value: str) -> bool:
    return len([x for x in value.split(";") if x.strip()]) == 2


def score_row(row: dict[str, str]) -> tuple:
    exact_name = int(
        norm_name(row.get("mibig_organism_name", ""))
        == norm_name(row.get("scientific_name", ""))
    )
    tax_match = int(
        row.get("mibig_ncbi_tax_id", "").strip()
        == row.get("tax_id", "").strip()
    )
    genomic = int(row.get("library_source", "").strip().upper() == "GENOMIC")
    random_sel = int(row.get("library_selection", "").strip().upper() == "RANDOM")
    two_fastqs = int(has_two_fastqs(row.get("fastq_ftp", "")))
    bases = parse_int(row.get("base_count", ""))
    reads = parse_int(row.get("read_count", ""))

    # Prefer complete exact metadata first; among those, prefer moderate data volume.
    # log-distance target is 1.5 Gbp, adequate for bacterial WGS while avoiding huge runs.
    target = 1_500_000_000
    distance = abs(math.log10(max(bases, 1)) - math.log10(target))

    return (
        exact_name,
        tax_match,
        genomic,
        random_sel,
        two_fastqs,
        -distance,
        bases,
        reads,
        row.get("run_accession", ""),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Select one deterministic primary ENA run per MIBiG locus from "
            "the broad paired Illumina WGS discovery table."
        )
    )
    parser.add_argument("--ena-tsv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    with args.ena_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []

    required = {
        "mibig_accession",
        "mibig_organism_name",
        "mibig_ncbi_tax_id",
        "run_accession",
        "scientific_name",
        "tax_id",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "base_count",
        "read_count",
        "fastq_ftp",
    }
    missing = required - set(fields)
    if missing:
        raise SystemExit("Missing ENA columns: " + ", ".join(sorted(missing)))

    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("library_strategy", "").strip().upper() != "WGS":
            continue
        if row.get("library_layout", "").strip().upper() != "PAIRED":
            continue
        if row.get("instrument_platform", "").strip().upper() != "ILLUMINA":
            continue
        if not has_two_fastqs(row.get("fastq_ftp", "")):
            continue
        groups[row["mibig_accession"].strip()].append(row)

    selected: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for accession in sorted(groups):
        candidates = groups[accession]
        exact = [
            row
            for row in candidates
            if norm_name(row.get("mibig_organism_name", ""))
            == norm_name(row.get("scientific_name", ""))
        ]
        pool = exact if exact else candidates
        chosen = max(pool, key=score_row)
        selected.append(chosen)
        summary_rows.append(
            {
                "mibig_accession": accession,
                "candidate_run_count": str(len(candidates)),
                "exact_normalized_name_run_count": str(len(exact)),
                "selection_basis": (
                    "EXACT_NORMALIZED_SCIENTIFIC_NAME"
                    if exact
                    else "TAXON_LEVEL_ONLY_REQUIRES_MANUAL_PROVENANCE_REVIEW"
                ),
                "selected_run_accession": chosen.get("run_accession", ""),
                "selected_scientific_name": chosen.get("scientific_name", ""),
                "selected_sample_accession": chosen.get("sample_accession", ""),
                "selected_base_count": chosen.get("base_count", ""),
                "selected_read_count": chosen.get("read_count", ""),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    summary_fields = [
        "mibig_accession",
        "candidate_run_count",
        "exact_normalized_name_run_count",
        "selection_basis",
        "selected_run_accession",
        "selected_scientific_name",
        "selected_sample_accession",
        "selected_base_count",
        "selected_read_count",
    ]
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    exact_count = sum(
        row["selection_basis"] == "EXACT_NORMALIZED_SCIENTIFIC_NAME"
        for row in summary_rows
    )
    manual_count = len(summary_rows) - exact_count
    print(f"Selected one primary run for {len(summary_rows)} MIBiG loci")
    print(f"Exact normalized-name loci: {exact_count}")
    print(f"Taxon-only/manual-review loci: {manual_count}")
    print(args.output)
    print(args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
