#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

DOWNLOAD_PROTOCOL_VERSION = "1.3-resumable-curl"
DEFAULT_ATTEMPTS_PER_TRANSPORT = 20


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


def md5sum(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256sum(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_curl(requested: str) -> Path:
    found = shutil.which(requested)
    if found is None:
        raise SystemExit(f"curl executable not found: {requested}")
    path = Path(found).resolve()
    if not path.is_file():
        raise SystemExit(f"Resolved curl executable is not a file: {path}")
    return path


def curl_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    first = (result.stdout or result.stderr).splitlines()
    return first[0].strip() if first else "UNKNOWN"


def curl_transfer(
    executable: Path,
    url: str,
    partial: Path,
    resume_from: int,
) -> subprocess.CompletedProcess[str]:
    command = [
        str(executable),
        "--location",
        "--fail",
        "--show-error",
        "--progress-bar",
        "--connect-timeout",
        "30",
        "--speed-limit",
        "1024",
        "--speed-time",
        "120",
        "--user-agent",
        "BGCPhaser-validation-download/1.3",
    ]
    if resume_from > 0:
        command.extend(["--continue-at", "-"])
    command.extend(["--output", str(partial), url])
    return subprocess.run(command, text=True)


def verified_existing(
    destination: Path,
    expected_md5: str,
    expected_size: int | None,
) -> dict[str, str] | None:
    if not destination.exists():
        return None
    observed_size = destination.stat().st_size
    if expected_size is not None and observed_size != expected_size:
        raise RuntimeError(
            f"Existing final file has wrong size and will not be overwritten: {destination} "
            f"(bytes={observed_size}, expected={expected_size})"
        )
    observed_md5 = md5sum(destination)
    if observed_md5.lower() != expected_md5.lower():
        raise RuntimeError(
            f"Existing final file has wrong MD5 and will not be overwritten: {destination} "
            f"(md5={observed_md5}, expected={expected_md5})"
        )
    return {
        "transport": "existing_verified",
        "requested_url": "",
        "observed_bytes": str(observed_size),
        "observed_md5": observed_md5,
    }


def validated_resumable_download(
    raw_url: str,
    destination: Path,
    expected_md5: str,
    expected_size: int | None,
    curl: Path,
    attempts_per_transport: int,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    diagnostics: list[dict[str, str]] = []

    existing = verified_existing(destination, expected_md5, expected_size)
    if existing is not None:
        return existing, diagnostics

    partial = destination.with_name(destination.name + ".part")
    if partial.exists() and expected_size is not None and partial.stat().st_size > expected_size:
        raise RuntimeError(
            f"Partial file is larger than ENA expected size; refusing to resume: {partial} "
            f"({partial.stat().st_size} > {expected_size})"
        )

    for transport, url in transport_urls(raw_url):
        for attempt in range(1, attempts_per_transport + 1):
            resume_from = partial.stat().st_size if partial.exists() else 0
            meta = {
                "transport": transport,
                "requested_url": url,
                "attempt": str(attempt),
                "resume_from_bytes": str(resume_from),
                "curl_exit_code": "",
                "observed_bytes": str(resume_from),
                "expected_bytes": "" if expected_size is None else str(expected_size),
                "observed_md5": "",
                "expected_md5": expected_md5,
                "result": "",
                "error": "",
            }

            print(
                f"[{destination.name}] {transport} attempt {attempt}/{attempts_per_transport}; "
                f"resume_from={resume_from} bytes",
                file=sys.stderr,
            )
            try:
                result = curl_transfer(curl, url, partial, resume_from)
            except Exception as exc:
                meta["result"] = "CURL_INVOCATION_ERROR"
                meta["error"] = repr(exc)
                diagnostics.append(meta)
                time.sleep(min(30, 3 * attempt))
                continue

            meta["curl_exit_code"] = str(result.returncode)
            observed_size = partial.stat().st_size if partial.exists() else 0
            meta["observed_bytes"] = str(observed_size)

            if expected_size is not None and observed_size > expected_size:
                meta["result"] = "OVERSIZE_PARTIAL"
                diagnostics.append(meta)
                raise RuntimeError(
                    f"Resumable download exceeded ENA expected size for {destination.name}: "
                    f"observed={observed_size}, expected={expected_size}. Partial retained at {partial}."
                )

            if result.returncode != 0:
                meta["result"] = "TRANSFER_INTERRUPTED_RETAINED"
                meta["error"] = f"curl_exit_code={result.returncode}"
                diagnostics.append(meta)
                # The partial file is intentionally retained. The next attempt
                # resumes from its current byte offset instead of starting over.
                time.sleep(min(30, 3 * attempt))
                continue

            if expected_size is not None and observed_size < expected_size:
                meta["result"] = "INCOMPLETE_TRANSFER_RETAINED"
                diagnostics.append(meta)
                time.sleep(min(30, 3 * attempt))
                continue

            observed_md5 = md5sum(partial)
            meta["observed_md5"] = observed_md5
            if observed_md5.lower() != expected_md5.lower():
                meta["result"] = "FULL_SIZE_MD5_MISMATCH"
                diagnostics.append(meta)
                raise RuntimeError(
                    f"Full-size file failed ENA MD5 for {destination.name}; partial retained for audit: "
                    f"observed={observed_md5}, expected={expected_md5}"
                )

            partial.replace(destination)
            meta["result"] = "VERIFIED"
            diagnostics.append(meta)
            return (
                {
                    "transport": f"curl_{transport}_resumable",
                    "requested_url": url,
                    "observed_bytes": str(observed_size),
                    "observed_md5": observed_md5,
                },
                diagnostics,
            )

    detail = "; ".join(
        f"{d['transport']}:{d['attempt']}:{d['result']}:"
        f"resume={d['resume_from_bytes']}:bytes={d['observed_bytes']}:curl={d['curl_exit_code'] or 'NA'}"
        for d in diagnostics[-10:]
    )
    raise RuntimeError(
        f"Resumable ENA download did not complete for {destination.name}. "
        f"Partial retained at {partial}. Last attempts: {detail}"
    )


def write_diagnostics(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "run_accession",
        "mate_index",
        "filename",
        "transport",
        "requested_url",
        "attempt",
        "resume_from_bytes",
        "curl_exit_code",
        "observed_bytes",
        "expected_bytes",
        "observed_md5",
        "expected_md5",
        "result",
        "error",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def preserve_prior_diagnostics(path: Path) -> Path | None:
    if not path.exists():
        return None
    candidate = path.with_name("download_diagnostics_pre_resume.tsv")
    if candidate.exists():
        index = 2
        while True:
            candidate = path.with_name(f"download_diagnostics_pre_resume_{index}.tsv")
            if not candidate.exists():
                break
            index += 1
    path.replace(candidate)
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Download exactly one paired-end ENA run using persistent curl resume, "
            "then verify ENA byte counts and MD5 checksums before finalizing files."
        )
    )
    parser.add_argument("--ena-tsv", required=True, type=Path)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--curl", default="curl")
    parser.add_argument(
        "--attempts-per-transport",
        type=int,
        default=DEFAULT_ATTEMPTS_PER_TRANSPORT,
    )
    args = parser.parse_args()

    if args.attempts_per_transport < 1:
        raise SystemExit("--attempts-per-transport must be >= 1")
    if not args.ena_tsv.is_file():
        raise SystemExit(f"ENA TSV absent: {args.ena_tsv}")

    curl = resolve_curl(args.curl)
    curl_ver = curl_version(curl)
    print(f"Downloader: {curl} ({curl_ver})", file=sys.stderr)

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

    preserved = preserve_prior_diagnostics(diagnostic_path)
    if preserved is not None:
        print(f"Preserved prior downloader diagnostics: {preserved}", file=sys.stderr)

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
            verified, attempts = validated_resumable_download(
                raw_url,
                destination,
                expected_md5,
                expected_size,
                curl,
                args.attempts_per_transport,
            )
            diagnostic_rows.extend(
                {
                    "run_accession": args.run,
                    "mate_index": str(index),
                    "filename": filename,
                    **attempt,
                }
                for attempt in attempts
            )
            print(
                f"Verified {filename}: bytes={verified['observed_bytes']} "
                f"md5={verified['observed_md5']}",
                file=sys.stderr,
            )
            manifest_rows.append(
                {
                    "run_accession": args.run,
                    "mate_index": str(index),
                    "filename": filename,
                    "transport": verified["transport"],
                    "ena_url": verified["requested_url"],
                    "ena_md5": expected_md5,
                    "observed_md5": verified["observed_md5"],
                    "expected_bytes": "" if expected_size is None else str(expected_size),
                    "observed_bytes": verified["observed_bytes"],
                    "sha256": sha256sum(destination),
                    "download_protocol": DOWNLOAD_PROTOCOL_VERSION,
                    "curl_executable": str(curl),
                    "curl_version": curl_ver,
                }
            )
            write_diagnostics(diagnostic_path, diagnostic_rows)
    except Exception as exc:
        if diagnostic_rows:
            write_diagnostics(diagnostic_path, diagnostic_rows)
        raise SystemExit(str(exc)) from exc

    write_diagnostics(diagnostic_path, diagnostic_rows)
    fields = [
        "run_accession",
        "mate_index",
        "filename",
        "transport",
        "ena_url",
        "ena_md5",
        "observed_md5",
        "expected_bytes",
        "observed_bytes",
        "sha256",
        "download_protocol",
        "curl_executable",
        "curl_version",
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
