#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

# Conservative discovery-stage flags. A match sends a run to manual review;
# it does not by itself exclude the run. Patterns target experimental contexts
# that can make a sequenced genome differ intentionally from the curated MIBiG
# reference strain or make the material unsuitable as ordinary isolate WGS.
DERIVATIVE_PATTERNS = {
    "MUTANT_OR_SUPPRESSOR": (
        r"\bmutants?\b|\bmutation(?:al|s)?\b|\bknock[- ]?out\b|"
        r"\bdeletion\b|\bdelta\b|\bsuppress(?:or|er)s?\b|\bsup\d*\b"
    ),
    "TRANSFORMANT_OR_ENGINEERED": (
        r"\btransformants?\b|\btransformation\b|\btransformed\b|"
        r"\bengineered\b|\bengineering\b|\bgenetically modified\b|"
        r"\brecombinant\b|\bsynthetic construct\b|\bminibacillus\b"
    ),
    "EVOLUTION_OR_PASSAGING": (
        r"\badaptive evolution\b|\blaboratory evolution\b|\bevolution experiment\b|"
        r"\bevolved\b|\bserial passage\b|\bpassaged\b|\bwhole population\b|"
        r"\bforeign DNA\b|\bgeneration\s*\d+\b|\bgen\s*\d+\b"
    ),
    "SELECTION_OR_RESISTANCE": (
        r"\bselected strain\b|\bselection experiment\b|\bresistan(?:t|ce)\b|"
        r"\bdrug[- ]selected\b|\bantibiotic[- ]selected\b"
    ),
    "IRRADIATION_OR_MUTAGENESIS": (
        r"\birradiat(?:ed|ion)\b|\bgamma ray\b|\bion beam\b|\bmutagenesis\b|"
        r"\bmock[- ]irradiated\b"
    ),
    "CRISPR_OR_REPORTER": (
        r"\bcrispr(?:i)?\b|\breporter strain\b|\breporter\b"
    ),
    "PHAGE_OR_PROPHAGE_EXPERIMENT": (
        r"\bprophages?\b|\btemperate phage\b|\bphage evolution\b|\blysogeny\b"
    ),
    "NONSTANDARD_CELL_DNA_SOURCE": (
        r"\bmembrane vesicles?\b|\bvesicle DNA\b"
    ),
    "CONTROL_OR_PARENT_CONTEXT": (
        r"\bancestor\b|\bparent strain\b|\bwild[- ]?type\b|\bWT background\b|"
        r"\bmock[- ]treated\b"
    ),
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
            "Conservatively flag ENA discovery runs whose study/sample/experiment "
            "metadata indicate experimental derivatives or nonstandard WGS contexts."
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
