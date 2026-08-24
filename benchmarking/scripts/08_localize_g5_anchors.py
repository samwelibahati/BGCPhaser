#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
from pathlib import Path

ANCHOR_LENGTH = 250
MIN_QUERY_COVERAGE = 0.95
MIN_IDENTITY = 0.98
MINIMAP2_PRESET = "asm5"


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_single_fasta(path: Path) -> tuple[str, str]:
    name: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError(f"Reference FASTA must contain exactly one sequence: {path}")
                name = line[1:].split()[0]
            else:
                if name is None:
                    raise ValueError(f"Sequence encountered before FASTA header: {path}")
                parts.append(line)
    if name is None or not parts:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return name, "".join(parts).upper()


def write_segment_fasta(gfa: Path, output: Path) -> tuple[int, int]:
    count = 0
    bases = 0
    with gfa.open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as target:
        for raw in source:
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed GFA S record: {raw.rstrip()}")
            segment_id, sequence = fields[1], fields[2]
            if sequence == "*":
                raise ValueError(f"GFA segment {segment_id} lacks sequence")
            target.write(f">{segment_id}\n{sequence}\n")
            count += 1
            bases += len(sequence)
    if count == 0:
        raise ValueError(f"No GFA segments found: {gfa}")
    return count, bases


def write_anchors(
    reference_name: str,
    sequence: str,
    start_1_inclusive: int,
    end_1_exclusive: int,
    output: Path,
) -> dict[str, tuple[int, int, str]]:
    if start_1_inclusive < 1:
        raise ValueError("start must be >= 1")
    if end_1_exclusive <= start_1_inclusive:
        raise ValueError("end must be greater than start")
    start0 = start_1_inclusive - 1
    end0 = end_1_exclusive - 1
    if end0 > len(sequence):
        raise ValueError(
            f"BGC end exceeds reference length: end0={end0}, reference={len(sequence)}"
        )
    if end0 - start0 < 2 * ANCHOR_LENGTH:
        raise ValueError("BGC interval is too short for two non-overlapping 250-nt anchors")

    anchors = {
        "left": (start0, start0 + ANCHOR_LENGTH, sequence[start0:start0 + ANCHOR_LENGTH]),
        "right": (end0 - ANCHOR_LENGTH, end0, sequence[end0 - ANCHOR_LENGTH:end0]),
    }
    with output.open("w", encoding="utf-8") as handle:
        for label, (left0, right0, anchor_sequence) in anchors.items():
            handle.write(
                f">{label}|reference={reference_name}|python0={left0}:{right0}|mibig1={left0 + 1}:{right0 + 1}\n"
                f"{anchor_sequence}\n"
            )
    return anchors


def resolve_executable(requested: str) -> Path:
    found = shutil.which(requested)
    if found is None:
        raise SystemExit(f"Executable not found: {requested}")
    path = Path(found).resolve()
    if not path.is_file():
        raise SystemExit(f"Resolved executable is not a file: {path}")
    return path


def minimap2_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"], check=True, capture_output=True, text=True
    )
    return result.stdout.strip() or result.stderr.strip()


def run_minimap2(executable: Path, target: Path, query: Path, paf: Path) -> None:
    command = [
        str(executable),
        "-x", MINIMAP2_PRESET,
        "--secondary=yes",
        "-N", "100",
        str(target),
        str(query),
    ]
    with paf.open("w", encoding="utf-8") as handle:
        subprocess.run(command, check=True, stdout=handle)


def parse_paf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            nmatch = int(fields[9])
            alnlen = int(fields[10])
            row = {
                "query": fields[0],
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": fields[4],
                "target": fields[5],
                "tlen": int(fields[6]),
                "tstart": int(fields[7]),
                "tend": int(fields[8]),
                "nmatch": nmatch,
                "alnlen": alnlen,
                "mapq": int(fields[11]),
                "query_coverage": (qend - qstart) / qlen if qlen else 0.0,
                "identity": nmatch / alnlen if alnlen else 0.0,
            }
            row["qualifies"] = (
                row["query_coverage"] >= MIN_QUERY_COVERAGE
                and row["identity"] >= MIN_IDENTITY
            )
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Prospective G5 BGC endpoint-anchor localization.")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--mibig-start", required=True, type=int, help="1-based inclusive")
    parser.add_argument("--mibig-end", required=True, type=int, help="1-based exclusive")
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--minimap2",
        default="minimap2",
        help="minimap2 executable name or explicit path; version 2.31 is required",
    )
    args = parser.parse_args()

    minimap2 = resolve_executable(args.minimap2)
    version = minimap2_version(minimap2)
    if not version.startswith("2.31"):
        raise SystemExit(f"Expected minimap2 2.31, observed {version}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    segments_fasta = args.output_dir / "gfa_segments.fasta"
    anchors_fasta = args.output_dir / "anchors_250nt.fasta"
    paf = args.output_dir / "anchor_alignments.paf"
    alignments_tsv = args.output_dir / "anchor_alignments.tsv"
    summary_tsv = args.output_dir / "g5_summary.tsv"

    reference_name, reference_sequence = read_single_fasta(args.reference_fasta)
    anchors = write_anchors(
        reference_name,
        reference_sequence,
        args.mibig_start,
        args.mibig_end,
        anchors_fasta,
    )
    segment_count, segment_bases = write_segment_fasta(args.gfa, segments_fasta)
    run_minimap2(minimap2, segments_fasta, anchors_fasta, paf)
    rows = parse_paf(paf)

    with alignments_tsv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "query", "target", "strand", "qlen", "qstart", "qend",
            "tlen", "tstart", "tend", "nmatch", "alnlen", "mapq",
            "query_coverage", "identity", "qualifies",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["query_coverage"] = f"{float(row['query_coverage']):.12f}"
            out["identity"] = f"{float(row['identity']):.12f}"
            writer.writerow(out)

    qualifying: dict[str, list[dict[str, object]]] = {"left": [], "right": []}
    for row in rows:
        label = str(row["query"]).split("|", 1)[0]
        if label in qualifying and bool(row["qualifies"]):
            qualifying[label].append(row)

    unique: dict[str, dict[str, object] | None] = {}
    for label in ("left", "right"):
        physical = {}
        for row in qualifying[label]:
            key = (str(row["target"]), str(row["strand"]), int(row["tstart"]), int(row["tend"]))
            physical[key] = row
        unique[label] = next(iter(physical.values())) if len(physical) == 1 else None

    status = "PASS"
    reason = ""
    if unique["left"] is None:
        status = "FAIL"
        reason = f"left_anchor_qualifying_localizations={len(qualifying['left'])}"
    elif unique["right"] is None:
        status = "FAIL"
        reason = f"right_anchor_qualifying_localizations={len(qualifying['right'])}"
    else:
        left = unique["left"]
        right = unique["right"]
        assert left is not None and right is not None
        if left["target"] == right["target"]:
            if left["strand"] != right["strand"]:
                status = "FAIL"
                reason = "same_segment_opposite_anchor_strands"
            elif left["strand"] == "+" and int(left["tstart"]) > int(right["tstart"]):
                status = "FAIL"
                reason = "same_segment_forward_coordinate_order_inconsistent"
            elif left["strand"] == "-" and int(left["tstart"]) < int(right["tstart"]):
                status = "FAIL"
                reason = "same_segment_reverse_coordinate_order_inconsistent"

    metrics: list[tuple[str, str]] = [
        ("benchmark_id", args.benchmark_id),
        ("reference_name", reference_name),
        ("reference_length", str(len(reference_sequence))),
        ("reference_sha256", sha256sum(args.reference_fasta)),
        ("gfa_sha256", sha256sum(args.gfa)),
        ("minimap2_executable", str(minimap2)),
        ("minimap2_sha256", sha256sum(minimap2)),
        ("minimap2_version", version),
        ("anchor_length", str(ANCHOR_LENGTH)),
        ("minimum_query_coverage", str(MIN_QUERY_COVERAGE)),
        ("minimum_identity", str(MIN_IDENTITY)),
        ("minimap2_preset", MINIMAP2_PRESET),
        ("mibig_start_1_inclusive", str(args.mibig_start)),
        ("mibig_end_1_exclusive", str(args.mibig_end)),
        ("gfa_segment_count", str(segment_count)),
        ("gfa_segment_bases", str(segment_bases)),
        ("left_anchor_reference_interval_python0", f"{anchors['left'][0]}:{anchors['left'][1]}"),
        ("right_anchor_reference_interval_python0", f"{anchors['right'][0]}:{anchors['right'][1]}"),
        ("left_qualifying_alignment_count", str(len(qualifying["left"]))),
        ("right_qualifying_alignment_count", str(len(qualifying["right"]))),
    ]
    for label in ("left", "right"):
        row = unique[label]
        if row is None:
            metrics.extend([
                (f"{label}_segment", "NA"),
                (f"{label}_orientation", "NA"),
                (f"{label}_target_interval", "NA"),
                (f"{label}_query_coverage", "NA"),
                (f"{label}_identity", "NA"),
            ])
        else:
            metrics.extend([
                (f"{label}_segment", str(row["target"])),
                (f"{label}_orientation", str(row["strand"])),
                (f"{label}_target_interval", f"{row['tstart']}:{row['tend']}"),
                (f"{label}_query_coverage", f"{float(row['query_coverage']):.12f}"),
                (f"{label}_identity", f"{float(row['identity']):.12f}"),
            ])
    metrics.extend([
        ("g5_status", status),
        ("g5_failure_reason", reason),
    ])

    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)

    for key, value in metrics:
        print(f"{key}\t{value}")
    print(f"alignments_tsv\t{alignments_tsv}")
    print(f"summary_tsv\t{summary_tsv}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
