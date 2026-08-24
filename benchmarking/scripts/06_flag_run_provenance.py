#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

DERIVATIVE_PATTERNS = {
    "MUTANT": r"\bmutant\b|\bmutation\b|\bknockout\b|\bknock-out\b|\bdeletion\b|\bdelta\b",
    "TRANSFORMANT": r"\btransformant\b|\btransformation\b|\btransformed\b",
    "ENGINEERED": r"\bengineered\b|\bgenetically modified\b|\brecombinant\b|\bsynthetic construct\b",
    "EVOLVED": r"\badaptive evolution\b|\blaboratory evolution\b|\bevolved strain\b|\bALE\b",
    "PASSAGED_OR_SELECTED": r"\bpassaged\b|\bserial passage\b|\bselected strain\b|\bselection experiment\b",
}

TEXT_FIELDS = [
    "study_title",
    "sample_title",
    "experiment_title",
    "library_name",
]


def norm(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip())


def provenance_flags(row: dict[str, str]) -> list[str]:
    text = " | ".join(norm(row.get(field, "")) for field in TEXT_FIELDS)
    flags: list[str] = []
    for label, pattern in DERIVATIVE_PATTERNS.items():
        if re.search(pattern, text, flags=re.IGNORECASE):
            flags.append(label)
    return flags


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Flag ENA discovery runs whose study/sample/experiment metadata "
            "suggests engineered, mutant, transformant, evolved or selected derivatives."
        )
    )
    parser.add_argument("--ena-tsv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.ena_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    required = {"mibig_accession", "run_accession", *TEXT_FIELDS}
    missing = required - set(fields)
    if missing:
        raise SystemExit("Missing provenance columns: " + ", ".join(sorted(missing)))

    out_fields = [*fields, "provenance_flag", "provenance_status"]
    flagged = 0
    clean = 0
    out_rows = []

    for row in rows:
        flags = provenance_flags(row)
        if flags:
            flagged += 1
            status = "REQUIRES_MANUAL_REVIEW"
        else:
            clean += 1
            status = "NO_DERIVATIVE_KEYWORD_DETECTED"
        out_rows.append(
            {
                **row,
                "provenance_flag": ";".join(flags),
                "provenance_status": status,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Rows: {len(out_rows)}")
    print(f"No derivative keyword detected: {clean}")
    print(f"Requires manual review: {flagged}")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
