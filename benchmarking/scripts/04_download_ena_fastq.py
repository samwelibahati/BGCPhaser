#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import time
import urllib.request
from pathlib import Path


def split_field(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def normalize_url(value: str) -> str:
    if value.startswith("ftp://"):
        return "https://" + value[len("ftp://"):]
    if value.startswith("http://") or value.startswith("https://"):
        return value
    return "https://" + value


def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, retries: int = 3) -> None:
    temporary = destination.with_suffix(destination.suffix + ".part")
    if temporary.exists():
        temporary.unlink()

    headers = {"User-Agent": "BGCPhaser-validation-download/1.0"}
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                with temporary.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            temporary.replace(destination)
            return
        except Exception as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(5 * attempt)

    raise RuntimeError(f"Download failed after {retries} attempts: {url}: {last_error}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Download exactly one paired-end ENA run from the TSV produced by "
            "02_query_ena_paired_wgs.py and verify ENA MD5 checksums."
        )
    )
    parser.add_argument("--ena-tsv", required=True, type=Path)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.ena_tsv.is_file():
        raise SystemExit(f"ENA TSV absent: {args.ena_tsv}")
    if args.output_dir.exists():
        raise SystemExit(f"Output directory already exists: {args.output_dir}")

    with args.ena_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader if row.get("run_accession", "").strip() == args.run]

    if len(rows) != 1:
        raise SystemExit(
            f"Expected exactly one ENA row for {args.run}; observed {len(rows)}"
        )

    row = rows[0]
    required_values = {
        "library_strategy": "WGS",
        "library_layout": "PAIRED",
        "instrument_platform": "ILLUMINA",
    }
    for field, expected in required_values.items():
        observed = row.get(field, "").strip().upper()
        if observed != expected:
            raise SystemExit(
                f"{args.run}: expected {field}={expected}, observed {observed or 'EMPTY'}"
            )

    urls = split_field(row.get("fastq_ftp", ""))
    md5s = split_field(row.get("fastq_md5", ""))
    sizes = split_field(row.get("fastq_bytes", ""))

    if len(urls) != 2:
        raise SystemExit(
            f"{args.run}: expected exactly two FASTQ files for paired data; observed {len(urls)}"
        )
    if len(md5s) != len(urls):
        raise SystemExit(
            f"{args.run}: FASTQ URL/MD5 cardinality mismatch: {len(urls)} vs {len(md5s)}"
        )
    if sizes and len(sizes) != len(urls):
        raise SystemExit(
            f"{args.run}: FASTQ URL/size cardinality mismatch: {len(urls)} vs {len(sizes)}"
        )

    args.output_dir.mkdir(parents=True)
    provenance = args.output_dir / "provenance.tsv"
    manifest_rows = []

    for index, (raw_url, expected_md5) in enumerate(zip(urls, md5s), start=1):
        url = normalize_url(raw_url)
        filename = Path(urllib.request.urlparse(url).path).name
        if not filename:
            raise SystemExit(f"Cannot derive filename from ENA URL: {url}")
        destination = args.output_dir / filename
        print(f"Downloading {args.run} mate {index}: {url}", file=sys.stderr)
        download(url, destination)

        observed_md5 = md5sum(destination)
        if observed_md5.lower() != expected_md5.lower():
            destination.unlink(missing_ok=True)
            raise SystemExit(
                f"MD5 mismatch for {filename}: expected {expected_md5}, observed {observed_md5}"
            )

        if sizes:
            expected_size = int(sizes[index - 1])
            observed_size = destination.stat().st_size
            if observed_size != expected_size:
                destination.unlink(missing_ok=True)
                raise SystemExit(
                    f"Size mismatch for {filename}: expected {expected_size}, observed {observed_size}"
                )
        else:
            observed_size = destination.stat().st_size

        manifest_rows.append(
            {
                "run_accession": args.run,
                "mate_index": str(index),
                "filename": filename,
                "ena_url": url,
                "ena_md5": expected_md5,
                "observed_md5": observed_md5,
                "bytes": str(observed_size),
                "sha256": sha256sum(destination),
            }
        )

    fields = [
        "run_accession",
        "mate_index",
        "filename",
        "ena_url",
        "ena_md5",
        "observed_md5",
        "bytes",
        "sha256",
    ]
    with provenance.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Verified {len(manifest_rows)} FASTQ files for {args.run}")
    print(provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
