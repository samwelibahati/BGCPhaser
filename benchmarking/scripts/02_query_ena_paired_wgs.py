#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ENA_ENDPOINT = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession",
    "study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "scientific_name",
    "tax_id",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]


def ena_search(tax_id: str) -> list[dict[str, str]]:
    query = f'tax_eq({tax_id}) AND library_strategy="WGS" AND library_layout="PAIRED" AND instrument_platform="ILLUMINA"'
    params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(FIELDS),
        "format": "json",
        "limit": "0",
    }
    url = ENA_ENDPOINT + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "BGCPhaser-validation-discovery/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = json.load(response)
    if not isinstance(payload, list):
        raise RuntimeError("Unexpected ENA response schema")
    rows = []
    for item in payload:
        if not isinstance(item, dict):
            continue
        rows.append({field: str(item.get(field, "")) for field in FIELDS})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Query ENA for paired-end Illumina WGS runs matching taxon IDs "
            "from a prospective MIBiG discovery pool."
        )
    )
    parser.add_argument("--mibig-pool", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--max-taxa",
        type=int,
        default=0,
        help="Optional development limit; 0 means all taxa.",
    )
    args = parser.parse_args()

    with args.mibig_pool.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"mibig_accession", "organism_name", "ncbi_tax_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit("Missing MIBiG pool columns: " + ", ".join(sorted(missing)))
        pool = list(reader)

    taxon_records: dict[str, list[dict[str, str]]] = {}
    for row in pool:
        tax_id = row["ncbi_tax_id"].strip()
        taxon_records.setdefault(tax_id, []).append(row)

    tax_ids = sorted(taxon_records)
    if args.max_taxa > 0:
        tax_ids = tax_ids[: args.max_taxa]

    output_rows: list[dict[str, str]] = []
    for index, tax_id in enumerate(tax_ids, start=1):
        print(f"[{index}/{len(tax_ids)}] ENA tax_id={tax_id}", file=sys.stderr)
        try:
            runs = ena_search(tax_id)
        except Exception as exc:
            print(f"ENA query failed for tax_id={tax_id}: {exc}", file=sys.stderr)
            continue

        for mibig in taxon_records[tax_id]:
            for run in runs:
                output_rows.append(
                    {
                        "mibig_accession": mibig["mibig_accession"],
                        "mibig_organism_name": mibig["organism_name"],
                        "mibig_ncbi_tax_id": tax_id,
                        **run,
                    }
                )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "mibig_accession",
        "mibig_organism_name",
        "mibig_ncbi_tax_id",
        *FIELDS,
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} MIBiG-run matches to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
