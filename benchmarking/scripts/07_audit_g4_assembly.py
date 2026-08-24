#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest(path: Path) -> tuple[int, int, list[str]]:
    total = 0
    passed = 0
    problems: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        match = re.match(r"^([0-9a-fA-F]{64})\s+(.+)$", line)
        if not match:
            problems.append(f"UNPARSEABLE:{line}")
            continue
        total += 1
        expected = match.group(1).lower()
        target = Path(match.group(2))
        if not target.is_absolute():
            target = (Path.cwd() / target).resolve()
        if not target.is_file():
            problems.append(f"MISSING:{target}")
            continue
        observed = sha256sum(target)
        if observed == expected:
            passed += 1
        else:
            problems.append(f"SHA256_MISMATCH:{target}:{expected}:{observed}")
    return total, passed, problems


def fasta_stats(path: Path) -> dict[str, int]:
    lengths: list[int] = []
    current = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    if not lengths:
        return {"count": 0, "total_bases": 0, "max_length": 0, "n50": 0}
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    cumulative = 0
    n50 = 0
    for length in ordered:
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return {
        "count": len(ordered),
        "total_bases": total,
        "max_length": max(ordered),
        "n50": n50,
    }


def gfa_stats(path: Path) -> dict[str, int]:
    segments = 0
    links = 0
    segment_bases = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("S\t"):
                segments += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[2] != "*":
                    segment_bases += len(fields[2])
            elif line.startswith("L\t"):
                links += 1
    return {"segments": segments, "links": links, "segment_bases": segment_bases}


def get_nested(data: dict, *keys: str):
    value = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return "NA"
        value = value[key]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit one completed BGCPhaser G4 assembly directory.")
    parser.add_argument("--g4-dir", required=True, type=Path)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--warnings-output", required=True, type=Path)
    args = parser.parse_args()

    g4 = args.g4_dir
    required = {
        "fastp_json": g4 / "fastp.json",
        "fastp_version": g4 / "logs/fastp.version.txt",
        "spades_version": g4 / "logs/spades.version.txt",
        "raw_manifest": g4 / "raw_checksums/raw_fastq.sha256",
        "processed_manifest": g4 / "reads/processed_fastq.sha256",
        "canonical_manifest": g4 / "spades/canonical_outputs.sha256",
        "gfa": g4 / "spades/assembly_graph_with_scaffolds.gfa",
        "contigs_paths": g4 / "spades/contigs.paths",
        "scaffolds_paths": g4 / "spades/scaffolds.paths",
        "contigs_fasta": g4 / "spades/contigs.fasta",
        "scaffolds_fasta": g4 / "spades/scaffolds.fasta",
        "spades_log": g4 / "spades/spades.log",
    }
    missing = [name for name, path in required.items() if not path.is_file()]
    if missing:
        raise SystemExit("Missing required G4 files: " + ", ".join(missing))

    audit: list[tuple[str, str]] = [("benchmark_id", args.benchmark_id)]

    fastp_version = required["fastp_version"].read_text(encoding="utf-8").strip()
    spades_version = required["spades_version"].read_text(encoding="utf-8").strip()
    audit.extend([
        ("fastp_version", fastp_version),
        ("spades_version", spades_version),
    ])

    fastp = json.loads(required["fastp_json"].read_text(encoding="utf-8"))
    before_reads = get_nested(fastp, "summary", "before_filtering", "total_reads")
    after_reads = get_nested(fastp, "summary", "after_filtering", "total_reads")
    before_bases = get_nested(fastp, "summary", "before_filtering", "total_bases")
    after_bases = get_nested(fastp, "summary", "after_filtering", "total_bases")
    audit.extend([
        ("fastp_before_total_reads", str(before_reads)),
        ("fastp_after_total_reads", str(after_reads)),
        ("fastp_before_total_bases", str(before_bases)),
        ("fastp_after_total_bases", str(after_bases)),
        ("fastp_after_q20_rate", str(get_nested(fastp, "summary", "after_filtering", "q20_rate"))),
        ("fastp_after_q30_rate", str(get_nested(fastp, "summary", "after_filtering", "q30_rate"))),
        ("fastp_after_gc_content", str(get_nested(fastp, "summary", "after_filtering", "gc_content"))),
        ("fastp_passed_filter_reads", str(get_nested(fastp, "filtering_result", "passed_filter_reads"))),
        ("fastp_low_quality_reads", str(get_nested(fastp, "filtering_result", "low_quality_reads"))),
        ("fastp_too_many_N_reads", str(get_nested(fastp, "filtering_result", "too_many_N_reads"))),
        ("fastp_too_short_reads", str(get_nested(fastp, "filtering_result", "too_short_reads"))),
    ])

    if isinstance(before_reads, int) and before_reads > 0 and isinstance(after_reads, int):
        audit.append(("fastp_read_retention_fraction", f"{after_reads / before_reads:.12f}"))
    if isinstance(before_bases, int) and before_bases > 0 and isinstance(after_bases, int):
        audit.append(("fastp_base_retention_fraction", f"{after_bases / before_bases:.12f}"))

    all_manifest_problems: list[str] = []
    for label, key in [
        ("raw", "raw_manifest"),
        ("processed", "processed_manifest"),
        ("canonical", "canonical_manifest"),
    ]:
        total, passed, problems = verify_manifest(required[key])
        audit.append((f"{label}_sha256_entries", str(total)))
        audit.append((f"{label}_sha256_verified", str(passed)))
        audit.append((f"{label}_sha256_status", "PASS" if total > 0 and total == passed and not problems else "FAIL"))
        all_manifest_problems.extend(f"{label}:{problem}" for problem in problems)

    contigs = fasta_stats(required["contigs_fasta"])
    scaffolds = fasta_stats(required["scaffolds_fasta"])
    graph = gfa_stats(required["gfa"])
    for prefix, stats in [("contigs", contigs), ("scaffolds", scaffolds)]:
        for key, value in stats.items():
            audit.append((f"{prefix}_{key}", str(value)))
    audit.extend([
        ("gfa_segments", str(graph["segments"])),
        ("gfa_links", str(graph["links"])),
        ("gfa_segment_bases", str(graph["segment_bases"])),
        ("gfa_bytes", str(required["gfa"].stat().st_size)),
        ("contigs_paths_bytes", str(required["contigs_paths"].stat().st_size)),
        ("scaffolds_paths_bytes", str(required["scaffolds_paths"].stat().st_size)),
    ])

    warning_lines: list[str] = []
    for line in required["spades_log"].read_text(encoding="utf-8", errors="replace").splitlines():
        lowered = line.lower()
        if "warning" in lowered or " warn" in lowered:
            warning_lines.append(line)
    args.warnings_output.parent.mkdir(parents=True, exist_ok=True)
    args.warnings_output.write_text(
        "\n".join(warning_lines) + ("\n" if warning_lines else ""),
        encoding="utf-8",
    )
    audit.append(("spades_warning_line_count", str(len(warning_lines))))
    audit.append(("sha256_manifest_problem_count", str(len(all_manifest_problems))))
    audit.append(("g4_integrity_status", "PASS" if not all_manifest_problems else "FAIL"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(audit)
        if all_manifest_problems:
            for problem in all_manifest_problems:
                writer.writerow(["sha256_problem", problem])

    for metric, value in audit:
        print(f"{metric}\t{value}")
    if all_manifest_problems:
        for problem in all_manifest_problems:
            print(f"sha256_problem\t{problem}")
    print(f"audit_output\t{args.output}")
    print(f"warnings_output\t{args.warnings_output}")
    return 0 if not all_manifest_problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
