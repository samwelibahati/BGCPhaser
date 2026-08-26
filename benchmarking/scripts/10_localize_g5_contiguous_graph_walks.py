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

import networkx as nx

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
ANALYSIS_ROLES = {
    "PRIMARY_VALIDATION",
    "DEVELOPMENT_SENTINEL_PROTOCOL_AMENDMENT_A1_1",
}


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_single_fasta(path: Path) -> tuple[str, str]:
    name: str | None = None
    sequence: list[str] = []
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
                    raise ValueError(f"Sequence before FASTA header: {path}")
                sequence.append(line)
    if name is None or not sequence:
        raise ValueError(f"No reference sequence found: {path}")
    return name, "".join(sequence).upper()


def write_anchors(
    reference_name: str,
    sequence: str,
    start_1_inclusive: int,
    end_1_exclusive: int,
    output: Path,
) -> dict[str, tuple[int, int, str]]:
    start0 = start_1_inclusive - 1
    end0 = end_1_exclusive - 1
    if start0 < 0 or end0 <= start0 or end0 > len(sequence):
        raise ValueError("Invalid MIBiG interval for reference sequence")
    if end0 - start0 < 2 * ANCHOR_LENGTH:
        raise ValueError("BGC interval shorter than two 250-nt anchors")
    anchors = {
        "left": (start0, start0 + ANCHOR_LENGTH, sequence[start0:start0 + ANCHOR_LENGTH]),
        "right": (end0 - ANCHOR_LENGTH, end0, sequence[end0 - ANCHOR_LENGTH:end0]),
    }
    with output.open("w", encoding="utf-8") as handle:
        for label, (left0, right0, anchor) in anchors.items():
            handle.write(
                f">{label}|reference={reference_name}|python0={left0}:{right0}"
                f"|mibig1={left0 + 1}:{right0 + 1}\n{anchor}\n"
            )
    return anchors


def resolve_executable(requested: str) -> Path:
    found = shutil.which(requested)
    if found is None:
        raise SystemExit(f"Executable not found: {requested}")
    resolved = Path(found).resolve()
    if not resolved.is_file():
        raise SystemExit(f"Resolved executable is not a file: {resolved}")
    return resolved


def minimap2_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() or result.stderr.strip()


def split_record_at_missing_edges(graph: nx.DiGraph, path):
    if not path:
        return [], []
    chunks = [[path[0]]]
    missing = []
    for left, right in zip(path, path[1:]):
        if graph.has_edge(left, right):
            chunks[-1].append(right)
        else:
            missing.append((left, right))
            chunks.append([right])
    return chunks, missing


def collect_contiguous_chunks(graph: nx.DiGraph, path_files: list[Path]):
    records_total = 0
    missing_transition_observations = 0
    unique_missing: set[tuple[tuple[str, str], tuple[str, str]]] = set()
    unique_chunks: dict[tuple[tuple[str, str], ...], list[tuple[str, str]]] = {}

    for path_file in path_files:
        records = parse_spades_path_records(path_file)
        records_total += len(records)
        for record in records:
            chunks, missing = split_record_at_missing_edges(graph, record)
            missing_transition_observations += len(missing)
            unique_missing.update(missing)
            for chunk in chunks:
                if chunk:
                    unique_chunks.setdefault(tuple(chunk), list(chunk))

    return (
        records_total,
        list(unique_chunks.values()),
        missing_transition_observations,
        unique_missing,
    )


def write_spelled_chunks(graph, chunks, fasta: Path, metadata: Path):
    chunk_by_id: dict[str, list[tuple[str, str]]] = {}
    sequence_by_id: dict[str, str] = {}
    with fasta.open("w", encoding="utf-8") as fasta_handle, metadata.open(
        "w", encoding="utf-8", newline=""
    ) as meta_handle:
        writer = csv.writer(meta_handle, delimiter="\t")
        writer.writerow(["chunk_id", "state_count", "spelled_length", "oriented_walk"])
        for index, chunk in enumerate(chunks, start=1):
            chunk_id = f"CW{index:06d}"
            spelled = spell_candidate(graph, chunk)
            chunk_by_id[chunk_id] = list(chunk)
            sequence_by_id[chunk_id] = spelled.sequence
            fasta_handle.write(f">{chunk_id}\n{spelled.sequence}\n")
            writer.writerow(
                [chunk_id, len(chunk), len(spelled.sequence), format_oriented_path(chunk)]
            )
    return chunk_by_id, sequence_by_id


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


def parse_paf(path: Path):
    rows = []
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


def reverse_oriented_path(path):
    return [(segment, flip_orientation(orientation)) for segment, orientation in reversed(path)]


def contribution_intervals(graph, path):
    spelled = spell_candidate(graph, path)
    intervals = []
    cursor = 0
    for index, state in enumerate(path):
        sequence = graph.nodes[state].get("sequence")
        if not isinstance(sequence, str):
            raise ValueError(f"Missing sequence for graph state {state}")
        overlap = 0 if index == 0 else spelled.overlaps[index - 1]
        global_start = cursor
        local_start = overlap
        global_end = cursor + len(sequence) - overlap
        intervals.append((global_start, global_end, local_start, len(sequence)))
        cursor = global_end
    if cursor != len(spelled.sequence):
        raise ValueError("Graph-walk coordinate construction failed")
    return spelled.sequence, intervals


def localize_alignment(graph, path, row):
    tlen = int(row["tlen"])
    strand = row["strand"]
    if strand == "+":
        oriented = list(path)
        aln_start = int(row["tstart"])
        aln_end = int(row["tend"])
    elif strand == "-":
        oriented = reverse_oriented_path(path)
        aln_start = tlen - int(row["tend"])
        aln_end = tlen - int(row["tstart"])
    else:
        raise ValueError(f"Unexpected PAF strand: {strand}")

    sequence, intervals = contribution_intervals(graph, oriented)
    if len(sequence) != tlen:
        raise ValueError("PAF target length differs from spelled graph-walk length")

    start_index = None
    end_index = None
    for index, (global_start, global_end, _local_start, _seg_len) in enumerate(intervals):
        if start_index is None and global_start <= aln_start < global_end:
            start_index = index
        if global_start <= aln_end - 1 < global_end:
            end_index = index
    if start_index is None or end_index is None:
        raise ValueError("Could not project anchor alignment onto graph states")

    start_global, _a, start_local, _b = intervals[start_index]
    end_global, _c, end_local, _d = intervals[end_index]
    start_offset = start_local + (aln_start - start_global)
    end_offset = end_local + (aln_end - end_global)
    subpath = oriented[start_index:end_index + 1]
    key = (tuple(subpath), int(start_offset), int(end_offset))
    return {
        "key": key,
        "subpath": subpath,
        "start_offset": int(start_offset),
        "end_offset": int(end_offset),
        "outer_start_state": subpath[0],
        "outer_end_state": subpath[-1],
    }


def state_text(state) -> str:
    return f"{state[0]}{state[1]}"


def primary_contribution(analysis_role: str, status: str) -> str:
    if analysis_role == "PRIMARY_VALIDATION":
        return "YES" if status == "PASS" else "NO_G5_FAILURE"
    return "NO_PROTOCOL_AMENDMENT_SENTINEL"


def main() -> int:
    parser = argparse.ArgumentParser(description="G5 A1.1 contiguous directed graph-walk anchor localization")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--contigs-paths", required=True, type=Path)
    parser.add_argument("--scaffolds-paths", required=True, type=Path)
    parser.add_argument("--mibig-start", required=True, type=int)
    parser.add_argument("--mibig-end", required=True, type=int)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument(
        "--analysis-role",
        required=True,
        choices=sorted(ANALYSIS_ROLES),
        help="Explicitly label this run as prospective primary validation or the frozen amendment sentinel.",
    )
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    for path in (args.reference_fasta, args.gfa, args.contigs_paths, args.scaffolds_paths):
        if not path.is_file():
            raise SystemExit(f"Required input missing: {path}")

    minimap2 = resolve_executable(args.minimap2)
    version = minimap2_version(minimap2)
    if not version.startswith("2.31"):
        raise SystemExit(f"Expected minimap2 2.31, observed {version}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    anchors_fasta = args.output_dir / "anchors_250nt.fasta"
    chunks_fasta = args.output_dir / "contiguous_graph_walks.fasta"
    chunks_tsv = args.output_dir / "contiguous_graph_walks.tsv"
    missing_tsv = args.output_dir / "split_transitions.tsv"
    paf = args.output_dir / "anchor_contiguous_walk_alignments.paf"
    alignments_tsv = args.output_dir / "anchor_contiguous_walk_alignments.tsv"
    localizations_tsv = args.output_dir / "g5_contiguous_walk_localizations.tsv"
    summary_tsv = args.output_dir / "g5_contiguous_walk_summary.tsv"

    reference_name, reference_sequence = read_single_fasta(args.reference_fasta)
    anchors = write_anchors(
        reference_name, reference_sequence, args.mibig_start, args.mibig_end, anchors_fasta
    )
    graph = read_gfa(args.gfa)

    records_total, chunks, missing_obs, unique_missing = collect_contiguous_chunks(
        graph, [args.contigs_paths, args.scaffolds_paths]
    )
    if not chunks:
        raise SystemExit("No contiguous graph-walk chunks were produced")

    with missing_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["left_state", "right_state"])
        for left, right in sorted(unique_missing):
            writer.writerow([state_text(left), state_text(right)])

    chunk_by_id, sequence_by_id = write_spelled_chunks(graph, chunks, chunks_fasta, chunks_tsv)
    run_minimap2(minimap2, chunks_fasta, anchors_fasta, paf)
    rows = parse_paf(paf)

    physical = {"left": {}, "right": {}}
    with alignments_tsv.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "query", "target_chunk_id", "strand", "qlen", "qstart", "qend",
            "tlen", "tstart", "tend", "nmatch", "alnlen", "mapq",
            "query_coverage", "identity", "qualifies", "physical_walk",
            "start_offset", "end_offset",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            target = str(row["target"])
            label = str(row["query"]).split("|", 1)[0]
            out = {
                "query": row["query"], "target_chunk_id": target, "strand": row["strand"],
                "qlen": row["qlen"], "qstart": row["qstart"], "qend": row["qend"],
                "tlen": row["tlen"], "tstart": row["tstart"], "tend": row["tend"],
                "nmatch": row["nmatch"], "alnlen": row["alnlen"], "mapq": row["mapq"],
                "query_coverage": f"{float(row['query_coverage']):.12f}",
                "identity": f"{float(row['identity']):.12f}",
                "qualifies": row["qualifies"], "physical_walk": "",
                "start_offset": "", "end_offset": "",
            }
            if bool(row["qualifies"]) and label in physical:
                loc = localize_alignment(graph, chunk_by_id[target], row)
                out["physical_walk"] = format_oriented_path(loc["subpath"])
                out["start_offset"] = loc["start_offset"]
                out["end_offset"] = loc["end_offset"]
                existing = physical[label].get(loc["key"])
                evidence = {
                    **loc,
                    "query_coverage": float(row["query_coverage"]),
                    "identity": float(row["identity"]),
                    "mapq": int(row["mapq"]),
                    "witness_chunks": {target},
                    "raw_alignment_count": 1,
                }
                if existing is None:
                    physical[label][loc["key"]] = evidence
                else:
                    existing["witness_chunks"].add(target)
                    existing["raw_alignment_count"] += 1
                    if (float(row["query_coverage"]), float(row["identity"]), int(row["mapq"])) > (
                        existing["query_coverage"], existing["identity"], existing["mapq"]
                    ):
                        existing["query_coverage"] = float(row["query_coverage"])
                        existing["identity"] = float(row["identity"])
                        existing["mapq"] = int(row["mapq"])
            writer.writerow(out)

    with localizations_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "anchor", "localization_id", "physical_walk", "state_count",
            "start_offset", "end_offset", "outer_start_state", "outer_end_state",
            "query_coverage", "identity", "mapq", "raw_alignment_count", "witness_chunks",
        ])
        for label in ("left", "right"):
            values = sorted(
                physical[label].values(),
                key=lambda x: (format_oriented_path(x["subpath"]), x["start_offset"], x["end_offset"]),
            )
            for index, item in enumerate(values, start=1):
                writer.writerow([
                    label, f"{label.upper()}_{index:03d}", format_oriented_path(item["subpath"]),
                    len(item["subpath"]), item["start_offset"], item["end_offset"],
                    state_text(item["outer_start_state"]), state_text(item["outer_end_state"]),
                    f"{item['query_coverage']:.12f}", f"{item['identity']:.12f}", item["mapq"],
                    item["raw_alignment_count"], ",".join(sorted(item["witness_chunks"])),
                ])

    left_values = list(physical["left"].values())
    right_values = list(physical["right"].values())
    status = "PASS"
    reason = ""
    if len(left_values) != 1:
        status = "FAIL"
        reason = f"left_anchor_physical_localizations={len(left_values)}"
    elif len(right_values) != 1:
        status = "FAIL"
        reason = f"right_anchor_physical_localizations={len(right_values)}"

    metrics = [
        ("benchmark_id", args.benchmark_id),
        ("analysis_role", args.analysis_role),
        ("reference_name", reference_name),
        ("reference_length", str(len(reference_sequence))),
        ("reference_sha256", sha256sum(args.reference_fasta)),
        ("gfa_sha256", sha256sum(args.gfa)),
        ("contigs_paths_sha256", sha256sum(args.contigs_paths)),
        ("scaffolds_paths_sha256", sha256sum(args.scaffolds_paths)),
        ("python_version", sys.version.split()[0]),
        ("networkx_version", nx.__version__),
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
        ("missing_edge_transition_observations", str(missing_obs)),
        ("unique_missing_edge_transitions", str(len(unique_missing))),
        ("unique_contiguous_graph_walk_chunks", str(len(chunks))),
        ("total_spelled_contiguous_walk_bases", str(sum(len(v) for v in sequence_by_id.values()))),
        ("left_physical_localization_count", str(len(left_values))),
        ("right_physical_localization_count", str(len(right_values))),
    ]

    for label, values in (("left", left_values), ("right", right_values)):
        if len(values) == 1:
            item = values[0]
            metrics.extend([
                (f"{label}_physical_walk", format_oriented_path(item["subpath"])),
                (f"{label}_start_offset", str(item["start_offset"])),
                (f"{label}_end_offset", str(item["end_offset"])),
                (f"{label}_outer_start_state", state_text(item["outer_start_state"])),
                (f"{label}_outer_end_state", state_text(item["outer_end_state"])),
                (f"{label}_query_coverage", f"{item['query_coverage']:.12f}"),
                (f"{label}_identity", f"{item['identity']:.12f}"),
            ])
        else:
            for suffix in (
                "physical_walk", "start_offset", "end_offset", "outer_start_state",
                "outer_end_state", "query_coverage", "identity",
            ):
                metrics.append((f"{label}_{suffix}", "NA"))

    if status == "PASS":
        left_item = left_values[0]
        right_item = right_values[0]
        metrics.extend([
            ("g6_outer_start_state", state_text(left_item["outer_start_state"])),
            ("g6_left_inner_state", state_text(left_item["outer_end_state"])),
            ("g6_right_inner_state", state_text(right_item["outer_start_state"])),
            ("g6_outer_end_state", state_text(right_item["outer_end_state"])),
        ])
    else:
        metrics.extend([
            ("g6_outer_start_state", "NA"), ("g6_left_inner_state", "NA"),
            ("g6_right_inner_state", "NA"), ("g6_outer_end_state", "NA"),
        ])

    metrics.extend([
        ("g5_contiguous_walk_status", status),
        ("g5_contiguous_walk_failure_reason", reason),
        ("primary_validation_contribution", primary_contribution(args.analysis_role, status)),
    ])

    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)

    for key, value in metrics:
        print(f"{key}\t{value}")
    print(f"split_transitions_tsv\t{missing_tsv}")
    print(f"alignments_tsv\t{alignments_tsv}")
    print(f"localizations_tsv\t{localizations_tsv}")
    print(f"summary_tsv\t{summary_tsv}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
