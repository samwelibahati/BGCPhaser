#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_field(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Lock one provenance-verified ENA run before validation download/assembly."
    )
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--mibig", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--expected-sample", required=True)
    parser.add_argument("--expected-study", required=True)
    parser.add_argument("--expected-reference", required=True)
    parser.add_argument("--expected-assembly", required=True)
    parser.add_argument(
        "--expected-library-selection",
        default="RANDOM",
        help=(
            "Exact ENA/SRA library_selection value expected for the prospectively "
            "selected run. Defaults to RANDOM to preserve prior lock behavior."
        ),
    )
    parser.add_argument(
        "--selection-basis",
        required=True,
        help="Pre-inspection provenance basis supporting this run selection.",
    )
    parser.add_argument(
        "--alternative-run",
        default="NA",
        help="Alternative same-source run considered before selection, if any.",
    )
    parser.add_argument(
        "--alternative-reason",
        default="NA",
        help="Prospectively recorded reason the alternative run was not selected.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.candidates.is_file():
        raise SystemExit(f"Candidate table absent: {args.candidates}")
    if not args.selection_basis.strip():
        raise SystemExit("--selection-basis must be non-empty")
    if not args.expected_library_selection.strip():
        raise SystemExit("--expected-library-selection must be non-empty")

    with args.candidates.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [
            row
            for row in reader
            if row.get("mibig_accession", "").strip() == args.mibig
            and row.get("run_accession", "").strip() == args.run
        ]

    if len(rows) != 1:
        raise SystemExit(
            f"Expected exactly one candidate row for {args.mibig}/{args.run}; observed {len(rows)}"
        )
    row = rows[0]

    checks = {
        "sample_accession": args.expected_sample,
        "study_accession": args.expected_study,
        "reference_accession": args.expected_reference,
        "library_strategy": "WGS",
        "library_layout": "PAIRED",
        "instrument_platform": "ILLUMINA",
        "library_source": "GENOMIC",
        "library_selection": args.expected_library_selection.strip(),
        "provenance_status": "NO_DERIVATIVE_KEYWORD_DETECTED",
    }
    failures: list[str] = []
    for field, expected in checks.items():
        observed = row.get(field, "").strip()
        if observed.upper() != expected.upper():
            failures.append(f"{field}: expected {expected!r}, observed {observed!r}")

    if row.get("provenance_flag", "").strip():
        failures.append(
            f"provenance_flag is non-empty: {row.get('provenance_flag', '').strip()!r}"
        )

    fastqs = split_field(row.get("fastq_ftp", ""))
    md5s = split_field(row.get("fastq_md5", ""))
    sizes = split_field(row.get("fastq_bytes", ""))
    if len(fastqs) != 2:
        failures.append(f"expected exactly two FASTQ URLs, observed {len(fastqs)}")
    if len(md5s) != 2:
        failures.append(f"expected exactly two FASTQ MD5 values, observed {len(md5s)}")
    if sizes and len(sizes) != 2:
        failures.append(f"expected zero or two FASTQ byte counts, observed {len(sizes)}")

    if failures:
        raise SystemExit("Validation-run lock failed:\n- " + "\n- ".join(failures))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ena_row_tsv = args.output_dir / "ena_run.tsv"
    selection_lock_tsv = args.output_dir / "selection_lock.tsv"
    if ena_row_tsv.exists() or selection_lock_tsv.exists():
        raise SystemExit(
            f"Lock outputs already exist in {args.output_dir}; refusing to overwrite"
        )

    # Preserve the exact discovery row so the downloader sees one and only one run record.
    with ena_row_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    lock_rows = [
        ("benchmark_id", args.benchmark_id),
        ("mibig_accession", args.mibig),
        ("mibig_organism_name", row.get("mibig_organism_name", "")),
        ("biosyn_class", row.get("biosyn_class", "")),
        ("reference_accession", args.expected_reference),
        ("assembly_accession", args.expected_assembly),
        ("sample_accession", args.expected_sample),
        ("secondary_sample_accession", row.get("secondary_sample_accession", "")),
        ("study_accession", args.expected_study),
        ("run_accession", args.run),
        ("experiment_accession", row.get("experiment_accession", "")),
        ("scientific_name", row.get("scientific_name", "")),
        ("study_title", row.get("study_title", "")),
        ("sample_title", row.get("sample_title", "")),
        ("experiment_title", row.get("experiment_title", "")),
        ("library_strategy", row.get("library_strategy", "")),
        ("library_source", row.get("library_source", "")),
        ("library_selection", row.get("library_selection", "")),
        ("expected_library_selection", args.expected_library_selection.strip()),
        ("library_layout", row.get("library_layout", "")),
        ("instrument_platform", row.get("instrument_platform", "")),
        ("instrument_model", row.get("instrument_model", "")),
        ("read_count", row.get("read_count", "")),
        ("base_count", row.get("base_count", "")),
        ("provenance_status", row.get("provenance_status", "")),
        ("provenance_flag", row.get("provenance_flag", "")),
        ("selection_basis", args.selection_basis.strip()),
        ("alternative_same_source_run", args.alternative_run.strip() or "NA"),
        ("alternative_reason_not_selected", args.alternative_reason.strip() or "NA"),
        ("score_inspected", "NO"),
        ("truth_rank_inspected", "NO"),
        ("chemistry_outcome_inspected", "NO"),
        ("reference_similarity_inspected", "NO"),
        ("candidate_table_sha256", sha256sum(args.candidates)),
        ("ena_run_tsv_sha256", sha256sum(ena_row_tsv)),
        ("lock_status", "PASS"),
    ]
    with selection_lock_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(lock_rows)

    for key, value in lock_rows:
        print(f"{key}\t{value}")
    print(f"ena_run_tsv\t{ena_row_tsv}")
    print(f"selection_lock_tsv\t{selection_lock_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
