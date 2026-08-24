#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bgcphaser.construction import spell_candidate
from bgcphaser.gfa import (
    flip_orientation,
    format_oriented_path,
    parse_spades_path_records,
    read_gfa,
)

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
                    raise ValueError(
                        f"Reference FASTA must contain exactly one sequence: {path}"
                    )
                name = line[1:].split()[0]
            else:
                if name is None:
                    raise ValueError(
                        f"Sequence encountered before FASTA header: {path}"
                    )
                parts.append(line)
    if name is None or not parts:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return name, "".join(parts).upper()


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
        "left": (
            start0,
            start0 + ANCHOR_LENGTH,
            sequence[start0 : start0 + ANCHOR_LENGTH],
        ),
        "right": (
            end0 - ANCHOR_LENGTH,
            end0,
            sequence[end0 - ANCHOR_LENGTH : end0],
        ),
    }
    with output.open("w", encoding="utf-8") as handle:
        for label, (left0, right0, anchor_sequence) in anchors.items():
            handle.write(
                f">{label}|reference={reference_name}|python0={left0}:{right0}"
                f"|mibig1={left0 + 1}:{right0 + 1}\n{anchor_sequence}\n"
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
        [str(executable), "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() or result.stderr.strip()


def read_unique_assembler_paths(path_files: list[Path]):
    records_total = 0
    unique: dict[tuple[tuple[str, str], ...], list[tuple[str, str]]] = {}
    for path_file in path_files:
        records = parse_spades_path_records(path_file)
        records_total += len(records)
        for path in records:
            key = tuple(path)
            unique.setdefault(key, list(path))
    return records_total, list(unique.values())


def write_spelled_paths(graph, paths, fasta: Path, metadata: Path):
    path_by_id: dict[str, list[tuple[str, str]]] = {}
    sequence_by_id: dict[str, str] = {}
    with fasta.open("w", encoding="utf-8") as fasta_handle, metadata.open(
        "w", encoding="utf-8", newline=""
    ) as meta_handle:
        writer = csv.writer(meta_handle, delimiter="\t")
        writer.writerow(["path_id", "state_count", "spelled_length", "oriented_path"])
        for index, path in enumerate(paths, start=1):
            path_id = f"AW{index:06d}"
            spelled = spell_candidate(graph, path)
            path_by_id[path_id] = list(path)
            sequence_by_id[path_id] = spelled.sequence
            fasta_handle.write(f">{path_id}\n{spelled.sequence}\n")
            writer.writerow(
                [path_id, len(path), len(spelled.sequence), format_oriented_path(path)]
            )
    return path_by_id, sequence_by_id


def run_minimap2(executable: Path, target: Path, query: Path, paf: Path) -> None:
    command = [
        str(executable),
        "-x",
        MINIMAP2_PRESET,
        "--secondary=yes",
        "-N",
        "100",
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
            row: dict[str, object] = {
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
                float(row["query_coverage"]) >= MIN_QUERY_COVERAGE
                and float(row["identity"]) >= MIN_IDENTITY
            )
            rows.append(row)
    return rows


def reverse_oriented_path(path: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return [
        (segment_id, flip_orientation(orientation))
        for segment_id, orientation in reversed(path)
    ]


def contribution_intervals(graph, path: list[tuple[str, str]]):
    spelled = spell_candidate(graph, path)
    intervals: list[tuple[int, int, int, int]] = []
    cursor = 0
    for index, state in enumerate(path):
        sequence = graph.nodes[state].get("sequence")
        if not isinstance(sequence, str):
            raise ValueError(f"Missing GFA sequence for state {state}")
        if index == 0:
            overlap = 0
        else:
            overlap = spelled.overlaps[index - 1]
        global_start = cursor
        local_start = overlap
        global_end = cursor + len(sequence) - overlap
        intervals.append((global_start, global_end, local_start, len(sequence)))
        cursor = global_end
    if cursor != len(spelled.sequence):
        raise ValueError("Contribution intervals do not span spelled path sequence")
    return spelled.sequence, intervals


def localize_alignment(graph, path: list[tuple[str, str]], row: dict[str, object]):
    tlen = int(row["tlen"])
    tstart = int(row["tstart"])
    tend = int(row["tend"])
    strand = str(row["strand"])

    if strand == "+":
        oriented = list(path)
        aln_start = tstart
        aln_end = tend
    elif strand == "-":
        oriented = reverse_oriented_path(path)
        aln_start = tlen - tend
        aln_end = tlen - tstart
    else:
        raise ValueError(f"Unexpected PAF strand: {strand}")

    spelled_sequence, intervals = contribution_intervals(graph, oriented)
    if len(spelled_sequence) != tlen:
        raise ValueError(
            f"PAF target length {tlen} differs from spelled path length {len(spelled_sequence)}"
        )
    if not (0 <= aln_start < aln_end <= tlen):
        raise ValueError(f"Alignment coordinates outside target: {aln_start}:{aln_end}/{tlen}")

    start_index = None
    end_index = None
    end_base = aln_end - 1
    for index, (global_start, global_end, _local_start, _segment_length) in enumerate(
        intervals
    ):
        if start_index is None and global_start <= aln_start < global_end:
            start_index = index
        if global_start <= end_base < global_end:
            end_index = index
    if start_index is None or end_index is None:
        raise ValueError("Could not project alignment boundaries onto graph-state contributions")
    if end_index < start_index:
        raise ValueError("Projected graph-state interval is reversed")

    start_global, _start_end, start_local_base, _start_len = intervals[start_index]
    end_global, _end_end, end_local_base, _end_len = intervals[end_index]
    start_offset = start_local_base + (aln_start - start_global)
    end_offset = end_local_base + (aln_end - end_global)
    subpath = oriented[start_index : end_index + 1]

    key = (
        tuple(subpath),
        int(start_offset),
        int(end_offset),
    )
    return {
        "key": key,
        "subpath": subpath,
        "start_offset": int(start_offset),
        "end_offset": int(end_offset),
        "outer_start_state": subpath[0],
        "outer_end_state": subpath[-1],
    }


def state_text(state: tuple[str, str]) -> str:
    return f"{state[0]}{state[1]}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Amended G5 graph-walk localization for predefined BGC endpoint anchors."
    )
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--contigs-paths", required=True, type=Path)
    parser.add_argument("--scaffolds-paths", required=True, type=Path)
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

    for required in (
        args.reference_fasta,
        args.gfa,
        args.contigs_paths,
        args.scaffolds_paths,
    ):
        if not required.is_file():
            raise SystemExit(f"Required input missing: {required}")

    minimap2 = resolve_executable(args.minimap2)
    version = minimap2_version(minimap2)
    if not version.startswith("2.31"):
        raise SystemExit(f"Expected minimap2 2.31, observed {version}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    anchors_fasta = args.output_dir / "anchors_250nt.fasta"
    paths_fasta = args.output_dir / "assembler_walks.fasta"
    paths_metadata = args.output_dir / "assembler_walks.tsv"
    paf = args.output_dir / "anchor_walk_alignments.paf"
    alignments_tsv = args.output_dir / "anchor_walk_alignments.tsv"
    localizations_tsv = args.output_dir / "g5_graph_walk_localizations.tsv"
    summary_tsv = args.output_dir / "g5_graph_walk_summary.tsv"

    reference_name, reference_sequence = read_single_fasta(args.reference_fasta)
    anchors = write_anchors(
        reference_name,
        reference_sequence,
        args.mibig_start,
        args.mibig_end,
        anchors_fasta,
    )

    graph = read_gfa(args.gfa)
    records_total, paths = read_unique_assembler_paths(
        [args.contigs_paths, args.scaffolds_paths]
    )
    if not paths:
        raise SystemExit("No SPAdes assembler paths were parsed")

    path_by_id, sequence_by_id = write_spelled_paths(
        graph, paths, paths_fasta, paths_metadata
    )
    run_minimap2(minimap2, paths_fasta, anchors_fasta, paf)
    rows = parse_paf(paf)

    physical_by_label: dict[
        str, dict[tuple[tuple[tuple[str, str], ...], int, int], dict[str, object]]
    ] = {"left": {}, "right": {}}

    with alignments_tsv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "query",
            "target_path_id",
            "strand",
            "qlen",
            "qstart",
            "qend",
            "tlen",
            "tstart",
            "tend",
            "nmatch",
            "alnlen",
            "mapq",
            "query_coverage",
            "identity",
            "qualifies",
            "physical_walk",
            "start_offset",
            "end_offset",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            target = str(row["target"])
            if target not in path_by_id:
                raise ValueError(f"PAF target absent from path metadata: {target}")
            out = {
                "query": row["query"],
                "target_path_id": target,
                "strand": row["strand"],
                "qlen": row["qlen"],
                "qstart": row["qstart"],
                "qend": row["qend"],
                "tlen": row["tlen"],
                "tstart": row["tstart"],
                "tend": row["tend"],
                "nmatch": row["nmatch"],
                "alnlen": row["alnlen"],
                "mapq": row["mapq"],
                "query_coverage": f"{float(row['query_coverage']):.12f}",
                "identity": f"{float(row['identity']):.12f}",
                "qualifies": row["qualifies"],
                "physical_walk": "",
                "start_offset": "",
                "end_offset": "",
            }
            label = str(row["query"]).split("|", 1)[0]
            if bool(row["qualifies"]) and label in physical_by_label:
                localization = localize_alignment(graph, path_by_id[target], row)
                out["physical_walk"] = format_oriented_path(localization["subpath"])
                out["start_offset"] = localization["start_offset"]
                out["end_offset"] = localization["end_offset"]
                key = localization["key"]
                current = physical_by_label[label].get(key)
                evidence = {
                    **localization,
                    "query_coverage": float(row["query_coverage"]),
                    "identity": float(row["identity"]),
                    "mapq": int(row["mapq"]),
                    "witness_path_ids": {target},
                    "raw_alignment_count": 1,
                }
                if current is None:
                    physical_by_label[label][key] = evidence
                else:
                    current["witness_path_ids"].add(target)
                    current["raw_alignment_count"] = int(current["raw_alignment_count"]) + 1
                    if (
                        float(row["query_coverage"]),
                        float(row["identity"]),
                        int(row["mapq"]),
                    ) > (
                        float(current["query_coverage"]),
                        float(current["identity"]),
                        int(current["mapq"]),
                    ):
                        current["query_coverage"] = float(row["query_coverage"])
                        current["identity"] = float(row["identity"])
                        current["mapq"] = int(row["mapq"])
            writer.writerow(out)

    with localizations_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "anchor",
                "localization_id",
                "physical_walk",
                "state_count",
                "start_offset",
                "end_offset",
                "outer_start_state",
                "outer_end_state",
                "query_coverage",
                "identity",
                "mapq",
                "raw_alignment_count",
                "witness_path_ids",
            ]
        )
        for label in ("left", "right"):
            items = sorted(
                physical_by_label[label].values(),
                key=lambda item: (
                    format_oriented_path(item["subpath"]),
                    int(item["start_offset"]),
                    int(item["end_offset"]),
                ),
            )
            for index, item in enumerate(items, start=1):
                writer.writerow(
                    [
                        label,
                        f"{label.upper()}_{index:03d}",
                        format_oriented_path(item["subpath"]),
                        len(item["subpath"]),
                        item["start_offset"],
                        item["end_offset"],
                        state_text(item["outer_start_state"]),
                        state_text(item["outer_end_state"]),
                        f"{float(item['query_coverage']):.12f}",
                        f"{float(item['identity']):.12f}",
                        item["mapq"],
                        item["raw_alignment_count"],
                        ",".join(sorted(item["witness_path_ids"])),
                    ]
                )

    left_values = list(physical_by_label["left"].values())
    right_values = list(physical_by_label["right"].values())
    status = "PASS"
    reason = ""
    if len(left_values) != 1:
        status = "FAIL"
        reason = f"left_anchor_physical_localizations={len(left_values)}"
    elif len(right_values) != 1:
        status = "FAIL"
        reason = f"right_anchor_physical_localizations={len(right_values)}"

    metrics: list[tuple[str, str]] = [
        ("benchmark_id", args.benchmark_id),
        ("analysis_role", "DEVELOPMENT_SENTINEL_PROTOCOL_AMENDMENT_A1"),
        ("reference_name", reference_name),
        ("reference_length", str(len(reference_sequence))),
        ("reference_sha256", sha256sum(args.reference_fasta)),
        ("gfa_sha256", sha256sum(args.gfa)),
        ("contigs_paths_sha256", sha256sum(args.contigs_paths)),
        ("scaffolds_paths_sha256", sha256sum(args.scaffolds_paths)),
        ("minimap2_executable", str(minimap2)),
        ("minimap2_sha256", sha256sum(minimap2)),
        ("minimap2_version", version),
        ("anchor_length", str(ANCHOR_LENGTH)),
        ("minimum_query_coverage", str(MIN_QUERY_COVERAGE)),
        ("minimum_identity", str(MIN_IDENTITY)),
        ("minimap2_preset", MINIMAP2_PRESET),
        ("mibig_start_1_inclusive", str(args.mibig_start)),
        ("mibig_end_1_exclusive", str(args.mibig_end)),
        ("left_anchor_reference_interval_python0", f"{anchors['left'][0]}:{anchors['left'][1]}"),
        ("right_anchor_reference_interval_python0", f"{anchors['right'][0]}:{anchors['right'][1]}"),
        ("assembler_path_records_total", str(records_total)),
        ("unique_oriented_assembler_paths", str(len(paths))),
        ("total_spelled_assembler_path_bases", str(sum(len(v) for v in sequence_by_id.values()))),
        ("left_physical_localization_count", str(len(left_values))),
        ("right_physical_localization_count", str(len(right_values))),
    ]

    for label, values in (("left", left_values), ("right", right_values)):
        if len(values) == 1:
            item = values[0]
            metrics.extend(
                [
                    (f"{label}_physical_walk", format_oriented_path(item["subpath"])),
                    (f"{label}_start_offset", str(item["start_offset"])),
                    (f"{label}_end_offset", str(item["end_offset"])),
                    (f"{label}_outer_start_state", state_text(item["outer_start_state"])),
                    (f"{label}_outer_end_state", state_text(item["outer_end_state"])),
                    (f"{label}_query_coverage", f"{float(item['query_coverage']):.12f}"),
                    (f"{label}_identity", f"{float(item['identity']):.12f}"),
                ]
            )
        else:
            metrics.extend(
                [
                    (f"{label}_physical_walk", "NA"),
                    (f"{label}_start_offset", "NA"),
                    (f"{label}_end_offset", "NA"),
                    (f"{label}_outer_start_state", "NA"),
                    (f"{label}_outer_end_state", "NA"),
                    (f"{label}_query_coverage", "NA"),
                    (f"{label}_identity", "NA"),
                ]
            )

    if status == "PASS":
        left_item = left_values[0]
        right_item = right_values[0]
        metrics.extend(
            [
                ("g6_outer_start_state", state_text(left_item["outer_start_state"])),
                ("g6_left_inner_state", state_text(left_item["outer_end_state"])),
                ("g6_right_inner_state", state_text(right_item["outer_start_state"])),
                ("g6_outer_end_state", state_text(right_item["outer_end_state"])),
            ]
        )
    else:
        metrics.extend(
            [
                ("g6_outer_start_state", "NA"),
                ("g6_left_inner_state", "NA"),
                ("g6_right_inner_state", "NA"),
                ("g6_outer_end_state", "NA"),
            ]
        )

    metrics.extend(
        [
            ("g5_graph_walk_status", status),
            ("g5_graph_walk_failure_reason", reason),
            ("primary_validation_contribution", "NO_PROTOCOL_AMENDMENT_SENTINEL"),
        ]
    )

    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)

    for key, value in metrics:
        print(f"{key}\t{value}")
    print(f"alignments_tsv\t{alignments_tsv}")
    print(f"localizations_tsv\t{localizations_tsv}")
    print(f"summary_tsv\t{summary_tsv}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
