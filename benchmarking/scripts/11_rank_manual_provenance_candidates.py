#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

DEFAULT_EXCLUDED_MIBIG = {"BGC0000309"}


def norm(value: str) -> str:
    text = value.strip().lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_int(value: str) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def two_fastqs(value: str) -> bool:
    return len([part for part in value.split(";") if part.strip()]) == 2


def metadata_text(row: dict[str, str]) -> str:
    fields = [
        "scientific_name",
        "study_title",
        "sample_title",
        "experiment_title",
        "library_name",
    ]
    return norm(" | ".join(row.get(field, "") for field in fields))


def tail_tokens(organism: str) -> list[str]:
    tokens = norm(organism).split()
    stop = {
        "strain",
        "str",
        "subsp",
        "ssp",
        "sp",
        "isolate",
        "bacterium",
    }
    # Treat the first two normalized tokens as the broad taxon label and use
    # the remainder only as a descriptive strain/isolate-specific signal.
    return [token for token in tokens[2:] if token not in stop and len(token) >= 2]


def tail_overlap(organism: str, metadata: str) -> tuple[int, int, str]:
    tokens = tail_tokens(organism)
    if not tokens:
        return 0, 0, ""
    metadata_tokens = set(metadata.split())
    matched = [token for token in tokens if token in metadata_tokens]
    return len(matched), len(tokens), ";".join(matched)


def ranking_key(row: dict[str, str]) -> tuple:
    clean = int(row["provenance_status"] == "NO_DERIVATIVE_KEYWORD_DETECTED")
    literal = int(row["mibig_organism_literal_in_metadata"] == "YES")
    exact = int(row["scientific_name_exact_normalized"] == "YES")
    tail_count = parse_int(row["organism_tail_token_match_count"])
    genomic = int(row["library_source"].strip().upper() == "GENOMIC")
    random_sel = int(row["library_selection"].strip().upper() == "RANDOM")
    paired_files = int(row["two_fastq_files"] == "YES")
    return (
        -clean,
        -literal,
        -exact,
        -tail_count,
        -genomic,
        -random_sel,
        -paired_files,
        row["mibig_accession"],
        row["run_accession"],
    )


def best_run_key(row: dict[str, str]) -> tuple:
    # Same provenance-first ordering as the global candidate table. Data
    # volume is reported for review but is deliberately not used to establish
    # isolate identity or automatic eligibility.
    return ranking_key(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank ENA/MIBiG pairs for manual provenance verification. This is "
            "a discovery aid only and never assigns validation eligibility."
        )
    )
    parser.add_argument("--mibig-pool", required=True, type=Path)
    parser.add_argument("--flagged-ena", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--one-per-locus", required=True, type=Path)
    parser.add_argument("--print-top", type=int, default=20)
    parser.add_argument(
        "--exclude-mibig",
        action="append",
        default=[],
        help="MIBiG accession to exclude; may be supplied repeatedly.",
    )
    args = parser.parse_args()

    excluded = set(DEFAULT_EXCLUDED_MIBIG)
    excluded.update(value.strip() for value in args.exclude_mibig if value.strip())

    with args.mibig_pool.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        pool_fields = reader.fieldnames or []
        required_pool = {
            "mibig_accession",
            "organism_name",
            "ncbi_tax_id",
            "biosyn_class",
            "reference_accession",
            "start_coord",
            "end_coord",
            "locus_span_bp",
            "locus_completeness",
            "total_module_count",
        }
        missing = required_pool - set(pool_fields)
        if missing:
            raise SystemExit("Missing MIBiG columns: " + ", ".join(sorted(missing)))
        pool = {row["mibig_accession"].strip(): row for row in reader}

    with args.flagged_ena.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        ena_fields = reader.fieldnames or []
        required_ena = {
            "mibig_accession",
            "mibig_organism_name",
            "run_accession",
            "study_accession",
            "study_title",
            "sample_accession",
            "secondary_sample_accession",
            "sample_title",
            "experiment_accession",
            "experiment_title",
            "scientific_name",
            "tax_id",
            "library_name",
            "library_strategy",
            "library_source",
            "library_selection",
            "library_layout",
            "instrument_platform",
            "instrument_model",
            "first_public",
            "read_count",
            "base_count",
            "fastq_ftp",
            "fastq_md5",
            "fastq_bytes",
            "provenance_flag",
            "provenance_status",
        }
        missing = required_ena - set(ena_fields)
        if missing:
            raise SystemExit("Missing ENA columns: " + ", ".join(sorted(missing)))
        ena_rows = list(reader)

    output_rows: list[dict[str, str]] = []
    for ena in ena_rows:
        accession = ena["mibig_accession"].strip()
        if accession in excluded or accession not in pool:
            continue
        mibig = pool[accession]

        if ena.get("library_strategy", "").strip().upper() != "WGS":
            continue
        if ena.get("library_layout", "").strip().upper() != "PAIRED":
            continue
        if ena.get("instrument_platform", "").strip().upper() != "ILLUMINA":
            continue
        if not two_fastqs(ena.get("fastq_ftp", "")):
            continue

        organism = mibig["organism_name"].strip()
        metadata = metadata_text(ena)
        organism_norm = norm(organism)
        scientific_norm = norm(ena.get("scientific_name", ""))
        matched_count, tail_count, matched_tokens = tail_overlap(organism, metadata)

        output_rows.append(
            {
                "review_status": "REQUIRES_MANUAL_PROVENANCE_VERIFICATION",
                "mibig_accession": accession,
                "biosyn_class": mibig.get("biosyn_class", ""),
                "mibig_organism_name": organism,
                "mibig_ncbi_tax_id": mibig.get("ncbi_tax_id", ""),
                "reference_accession": mibig.get("reference_accession", ""),
                "start_coord": mibig.get("start_coord", ""),
                "end_coord": mibig.get("end_coord", ""),
                "locus_span_bp": mibig.get("locus_span_bp", ""),
                "locus_completeness": mibig.get("locus_completeness", ""),
                "total_module_count": mibig.get("total_module_count", ""),
                "run_accession": ena.get("run_accession", ""),
                "study_accession": ena.get("study_accession", ""),
                "study_title": ena.get("study_title", ""),
                "sample_accession": ena.get("sample_accession", ""),
                "secondary_sample_accession": ena.get("secondary_sample_accession", ""),
                "sample_title": ena.get("sample_title", ""),
                "experiment_accession": ena.get("experiment_accession", ""),
                "experiment_title": ena.get("experiment_title", ""),
                "scientific_name": ena.get("scientific_name", ""),
                "tax_id": ena.get("tax_id", ""),
                "library_name": ena.get("library_name", ""),
                "library_strategy": ena.get("library_strategy", ""),
                "library_source": ena.get("library_source", ""),
                "library_selection": ena.get("library_selection", ""),
                "library_layout": ena.get("library_layout", ""),
                "instrument_platform": ena.get("instrument_platform", ""),
                "instrument_model": ena.get("instrument_model", ""),
                "first_public": ena.get("first_public", ""),
                "read_count": ena.get("read_count", ""),
                "base_count": ena.get("base_count", ""),
                "two_fastq_files": "YES",
                "scientific_name_exact_normalized": (
                    "YES" if organism_norm == scientific_norm else "NO"
                ),
                "mibig_organism_literal_in_metadata": (
                    "YES" if organism_norm and organism_norm in metadata else "NO"
                ),
                "organism_tail_token_match_count": str(matched_count),
                "organism_tail_token_count": str(tail_count),
                "organism_tail_tokens_matched": matched_tokens,
                "provenance_status": ena.get("provenance_status", ""),
                "provenance_flag": ena.get("provenance_flag", ""),
                "fastq_ftp": ena.get("fastq_ftp", ""),
                "fastq_md5": ena.get("fastq_md5", ""),
                "fastq_bytes": ena.get("fastq_bytes", ""),
            }
        )

    output_rows.sort(key=ranking_key)

    fields = [
        "review_status",
        "mibig_accession",
        "biosyn_class",
        "mibig_organism_name",
        "mibig_ncbi_tax_id",
        "reference_accession",
        "start_coord",
        "end_coord",
        "locus_span_bp",
        "locus_completeness",
        "total_module_count",
        "run_accession",
        "study_accession",
        "study_title",
        "sample_accession",
        "secondary_sample_accession",
        "sample_title",
        "experiment_accession",
        "experiment_title",
        "scientific_name",
        "tax_id",
        "library_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "first_public",
        "read_count",
        "base_count",
        "two_fastq_files",
        "scientific_name_exact_normalized",
        "mibig_organism_literal_in_metadata",
        "organism_tail_token_match_count",
        "organism_tail_token_count",
        "organism_tail_tokens_matched",
        "provenance_status",
        "provenance_flag",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in output_rows:
        groups[row["mibig_accession"]].append(row)

    one_per_locus = [min(rows, key=best_run_key) for rows in groups.values()]
    one_per_locus.sort(key=ranking_key)
    with args.one_per_locus.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(one_per_locus)

    print(f"Excluded MIBiG accessions: {','.join(sorted(excluded))}")
    print(f"Candidate run-locus pairs: {len(output_rows)}")
    print(f"Candidate loci represented: {len(one_per_locus)}")
    print("No row is automatically accepted; all require manual provenance verification.")
    print()

    top_n = max(0, args.print_top)
    print("rank\tmibig_accession\tmibig_organism_name\treference_accession\trun_accession\tsample_accession\tstudy_accession\tprovenance_status\texact_name\tliteral_organism\ttail_match\tbase_count\tstudy_title\tsample_title\texperiment_title")
    for index, row in enumerate(one_per_locus[:top_n], start=1):
        tail = f"{row['organism_tail_token_match_count']}/{row['organism_tail_token_count']}"
        values = [
            str(index),
            row["mibig_accession"],
            row["mibig_organism_name"],
            row["reference_accession"],
            row["run_accession"],
            row["sample_accession"],
            row["study_accession"],
            row["provenance_status"],
            row["scientific_name_exact_normalized"],
            row["mibig_organism_literal_in_metadata"],
            tail,
            row["base_count"],
            row["study_title"],
            row["sample_title"],
            row["experiment_title"],
        ]
        print("\t".join(value.replace("\t", " ").replace("\n", " ") for value in values))

    print()
    print(args.output)
    print(args.one_per_locus)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
