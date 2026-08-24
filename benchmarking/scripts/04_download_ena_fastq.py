#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


def split_field(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def transport_urls(value: str) -> list[tuple[str, str]]:
    raw = value.strip()
    if raw.startswith("https://"):
        hostpath = raw[len("https://"):]
    elif raw.startswith("http://"):
        hostpath = raw[len("http://"):]
    elif raw.startswith("ftp://"):
        hostpath = raw[len("ftp://"):]
    else:
        hostpath = raw
    return [("https", "https://" + hostpath), ("ftp", "ftp://" + hostpath)]


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


def download_once(url: str, temporary: Path) -> dict[str, str]:
    headers = {"User-Agent": "BGCPhaser-validation-download/1.2"}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=180) as response:
        final_url = response.geturl()
        status = getattr(response, "status", "")
        content_type = response.headers.get("Content-Type", "")
        content_length = response.headers.get("Content-Length", "")
        with temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    return {
        "final_url": str(final_url),
        "http_status": str(status),
        "content_type": str(content_type),
        "content_length": str(content_length),
    }


def validated_download(
    raw_url: str,
    destination: Path,
    expected_md5: str,
    expected_size: int | None,
    retries: int = 5,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    diagnostics: list[dict[str, str]] = []

    if destination.exists():
        observed_size = destination.stat().st_size
        observed_md5 = md5sum(destination)
        size_ok = expected_size is None or observed_size == expected_size
        md5_ok = observed_md5.lower() == expected_md5.lower()
        if size_ok and md5_ok:
            return (
                {
                    "transport": "existing_verified",
                    "requested_url": "",
                    "final_url": "",
                    "http_status": "",
                    "content_type": "",
                    "content_length": "",
                    "observed_bytes": str(observed_size),
                    "observed_md5": observed_md5,
                },
                diagnostics,
            )
        raise RuntimeError(
            f"Existing file is invalid and will not be overwritten: {destination} "
            f"(bytes={observed_size}, md5={observed_md5})"
        )

    for transport, url in transport_urls(raw_url):
        for attempt in range(1, retries + 1):
            temporary = destination.with_name(
                destination.name + f".{transport}.attempt{attempt}.part"
            )
            temporary.unlink(missing_ok=True)
            meta = {
                "transport": transport,
                "requested_url": url,
                "attempt": str(attempt),
                "final_url": "",
                "http_status": "",
                "content_type": "",
                "content_length": "",
                "observed_bytes": "",
                "expected_bytes": "" if expected_size is None else str(expected_size),
                "observed_md5": "",
                "expected_md5": expected_md5,
                "result": "",
                "error": "",
            }
            try:
                response_meta = download_once(url, temporary)
                meta.update(response_meta)
                observed_size = temporary.stat().st_size
                observed_md5 = md5sum(temporary)
                meta["observed_bytes"] = str(observed_size)
                meta["observed_md5"] = observed_md5

                if expected_size is not None and observed_size != expected_size:
                    meta["result"] = "SIZE_MISMATCH"
                    diagnostics.append(meta)
                    temporary.unlink(missing_ok=True)
                    if attempt < retries:
                        time.sleep(5 * attempt)
                    continue

                if observed_md5.lower() != expected_md5.lower():
                    meta["result"] = "MD5_MISMATCH"
                    diagnostics.append(meta)
                    temporary.unlink(missing_ok=True)
                    if attempt < retries:
                        time.sleep(5 * attempt)
                    continue

                temporary.replace(destination)
                meta["result"] = "VERIFIED"
                diagnostics.append(meta)
                return (
                    {
                        "transport": transport,
                        "requested_url": url,
                        "final_url": meta["final_url"],
                        "http_status": meta["http_status"],
                        "content_type": meta["content_type"],
                        "content_length": meta["content_length"],
                        "observed_bytes": str(observed_size),
                        "observed_md5": observed_md5,
                    },
                    diagnostics,
                )
            except Exception as exc:
                meta["result"] = "TRANSFER_ERROR"
                meta["error"] = repr(exc)
                diagnostics.append(meta)
                temporary.unlink(missing_ok=True)
                if attempt < retries:
                    time.sleep(5 * attempt)

    detail = "; ".join(
        f"{d['transport']}:{d['attempt']}:{d['result']}:"
        f"bytes={d['observed_bytes'] or 'NA'}:md5={d['observed_md5'] or 'NA'}"
        for d in diagnostics
    )
    raise RuntimeError(
        f"No ENA transport produced a verified file for {destination.name}: {detail}"
    )


def write_diagnostics(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "run_accession", "mate_index", "filename", "transport", "requested_url",
        "attempt", "final_url", "http_status", "content_type", "content_length",
        "observed_bytes", "expected_bytes", "observed_md5", "expected_md5",
        "result", "error",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Download exactly one paired-end ENA run and verify ENA byte counts "
            "and MD5 checksums. HTTPS is attempted first, followed by FTP."
        )
    )
    parser.add_argument("--ena-tsv", required=True, type=Path)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.ena_tsv.is_file():
        raise SystemExit(f"ENA TSV absent: {args.ena_tsv}")

    with args.ena_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader if row.get("run_accession", "").strip() == args.run]

    if len(rows) != 1:
        raise SystemExit(f"Expected exactly one ENA row for {args.run}; observed {len(rows)}")

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
        raise SystemExit(f"{args.run}: expected exactly two FASTQ files; observed {len(urls)}")
    if len(md5s) != len(urls):
        raise SystemExit(f"{args.run}: FASTQ URL/MD5 cardinality mismatch")
    if sizes and len(sizes) != len(urls):
        raise SystemExit(f"{args.run}: FASTQ URL/size cardinality mismatch")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    provenance = args.output_dir / "provenance.tsv"
    diagnostic_path = args.output_dir / "download_diagnostics.tsv"
    if provenance.exists():
        raise SystemExit(f"Verified provenance already exists: {provenance}")

    manifest_rows: list[dict[str, str]] = []
    diagnostic_rows: list[dict[str, str]] = []

    try:
        for index, (raw_url, expected_md5) in enumerate(zip(urls, md5s), start=1):
            parsed = urllib.parse.urlparse(
                raw_url if "://" in raw_url else "https://" + raw_url
            )
            filename = Path(parsed.path).name
            if not filename:
                raise RuntimeError(f"Cannot derive filename from ENA URL: {raw_url}")
            destination = args.output_dir / filename
            expected_size = int(sizes[index - 1]) if sizes else None

            print(
                f"Downloading {args.run} mate {index}; expected bytes={expected_size or 'NA'}; "
                f"expected md5={expected_md5}",
                file=sys.stderr,
            )
            verified, attempts = validated_download(
                raw_url, destination, expected_md5, expected_size
            )
            for attempt in attempts:
                diagnostic_rows.append(
                    {
                        "run_accession": args.run,
                        "mate_index": str(index),
                        "filename": filename,
                        **attempt,
                    }
                )
            print(
                f"Verified {filename} via {verified['transport']}: "
                f"bytes={verified['observed_bytes']} md5={verified['observed_md5']}",
                file=sys.stderr,
            )
            manifest_rows.append(
                {
                    "run_accession": args.run,
                    "mate_index": str(index),
                    "filename": filename,
                    "transport": verified["transport"],
                    "ena_url": verified["requested_url"],
                    "final_url": verified["final_url"],
                    "ena_md5": expected_md5,
                    "observed_md5": verified["observed_md5"],
                    "expected_bytes": "" if expected_size is None else str(expected_size),
                    "observed_bytes": verified["observed_bytes"],
                    "sha256": sha256sum(destination),
                }
            )
    except Exception as exc:
        if diagnostic_rows:
            write_diagnostics(diagnostic_path, diagnostic_rows)
        raise SystemExit(str(exc)) from exc

    write_diagnostics(diagnostic_path, diagnostic_rows)
    fields = [
        "run_accession", "mate_index", "filename", "transport", "ena_url",
        "final_url", "ena_md5", "observed_md5", "expected_bytes",
        "observed_bytes", "sha256",
    ]
    with provenance.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Verified {len(manifest_rows)} FASTQ files for {args.run}")
    print(provenance)
    print(diagnostic_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
