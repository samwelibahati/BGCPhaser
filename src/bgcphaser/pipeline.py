from __future__ import annotations

from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
from typing import Mapping, Sequence

import networkx as nx

from .chemistry import ChemistryFeatures, chemistry_from_antismash_json
from .construction import (
    TransitionClasses,
    candidate_transition_support,
    construct_transition_classes,
    parse_oriented_path,
    spell_candidate,
    write_fasta,
)
from .evidence import (
    grouped_sam_records,
    l_support_decision,
    r_support_decision,
    support_counts,
)
from .gfa import (
    OrientedNode,
    flip_orientation,
    format_oriented_path,
    parse_spades_path_records,
    read_gfa,
)
from .scoring import format_q12, rank_candidates, score_candidate
from .sequence import (
    AssemblerDiscordance,
    CoverageConsistency,
    SequenceComposite,
    assembler_discordance,
    coverage_consistency,
    sequence_composite,
)
from .tools import run_antismash, run_minimap2


R_CONTEXT_BP = 30


@dataclass(frozen=True)
class AnalysisRun:
    output_dir: Path
    ranked_path: Path
    candidate_count: int
    ranked_count: int


def _require_file(path: str | Path, *, label: str) -> Path:
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} absent: {resolved}")
    return resolved


def _read_candidates(path: str | Path) -> dict[str, list[OrientedNode]]:
    candidate_path = _require_file(path, label="candidate table")
    with candidate_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"candidate_id", "oriented_walk"}
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(
                "Candidate table missing columns: " + ", ".join(missing)
            )
        candidates: dict[str, list[OrientedNode]] = {}
        for row_number, row in enumerate(reader, start=2):
            candidate_id = row["candidate_id"].strip()
            if not candidate_id:
                raise ValueError(
                    f"Candidate row {row_number}: candidate_id is empty"
                )
            if candidate_id in candidates:
                raise ValueError(f"Duplicate candidate_id: {candidate_id}")
            try:
                candidates[candidate_id] = parse_oriented_path(
                    row["oriented_walk"]
                )
            except ValueError as exc:
                raise ValueError(
                    f"Candidate row {row_number}: {exc}"
                ) from exc
    if not candidates:
        raise ValueError("Candidate table contains no candidates")
    return candidates


def _write_tsv(path: Path, fields: Sequence[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _depth_by_segment(
    graph: nx.DiGraph,
    depth_tag: str | None,
) -> dict[str, float]:
    if depth_tag is None:
        return {}
    depth_tag = depth_tag.strip()
    if not depth_tag:
        raise ValueError("depth_tag cannot be empty")

    result: dict[str, float] = {}
    for state, attributes in graph.nodes(data=True):
        segment_id, orientation = state
        if orientation != "+":
            continue
        values = []
        for raw_tag in attributes.get("raw_tags", ()):
            pieces = raw_tag.split(":", 2)
            if len(pieces) == 3 and pieces[0] == depth_tag:
                values.append(pieces[2])
        if len(values) > 1:
            raise ValueError(
                f"Segment {segment_id} has duplicate {depth_tag} tags"
            )
        if not values:
            continue
        try:
            value = float(values[0])
        except ValueError as exc:
            raise ValueError(
                f"Segment {segment_id} has invalid {depth_tag} value: {values[0]!r}"
            ) from exc
        if not math.isfinite(value):
            raise ValueError(
                f"Segment {segment_id} has non-finite {depth_tag} value"
            )
        result[segment_id] = value

    if not result:
        raise ValueError(
            f"No usable GFA depth tags named {depth_tag!r} were found"
        )
    return result


def _reverse_state(state: OrientedNode) -> OrientedNode:
    return state[0], flip_orientation(state[1])


def _assembler_adjacency(
    path_files: Sequence[str | Path],
) -> set[tuple[OrientedNode, OrientedNode]] | None:
    if not path_files:
        return None
    adjacency: set[tuple[OrientedNode, OrientedNode]] = set()
    for raw_path in path_files:
        path = _require_file(raw_path, label="assembler path file")
        for record in parse_spades_path_records(path):
            for left, right in zip(record, record[1:]):
                adjacency.add((left, right))
                adjacency.add((_reverse_state(right), _reverse_state(left)))
    return adjacency


def _r_support_from_sam(
    sam_path: Path,
    target_lengths: Mapping[str, int],
) -> Counter[str]:
    decisions = []
    with sam_path.open("r", encoding="utf-8") as handle:
        for _query_name, records in grouped_sam_records(handle):
            first = []
            second = []
            unpaired = []
            for record in records:
                if record.first and record.second:
                    raise ValueError(
                        "SAM record is marked as both first and second mate"
                    )
                if record.first:
                    first.append(record)
                elif record.second:
                    second.append(record)
                else:
                    unpaired.append(record)
            if unpaired and (first or second):
                raise ValueError(
                    "SAM query mixes paired and unpaired record flags"
                )
            for read_end_records in (first, second, unpaired):
                if read_end_records:
                    decisions.append(
                        r_support_decision(read_end_records, target_lengths)
                    )
    return support_counts(decisions)


def _l_support_from_sam(
    sam_path: Path,
    *,
    junction_by_class: Mapping[str, int],
    insert_min: int,
    insert_max: int,
) -> Counter[str]:
    decisions = []
    saw_paired_query = False
    with sam_path.open("r", encoding="utf-8") as handle:
        for _query_name, records in grouped_sam_records(handle):
            has_first = any(record.first for record in records)
            has_second = any(record.second for record in records)
            if has_first and has_second:
                saw_paired_query = True
            decision, _geometry = l_support_decision(
                records,
                junction_by_class=junction_by_class,
                insert_min=insert_min,
                insert_max=insert_max,
            )
            decisions.append(decision)
    if junction_by_class and not saw_paired_query:
        raise ValueError(
            "L evidence requires paired-end SAM first/second mate flags"
        )
    return support_counts(decisions)


def _candidate_sequence_classes(
    graph: nx.DiGraph,
    candidates: Mapping[str, Sequence[OrientedNode]],
):
    representative_by_digest: dict[str, str] = {}
    sequence_by_digest: dict[str, str] = {}
    digest_by_candidate: dict[str, str] = {}

    for candidate_id in sorted(candidates):
        sequence = spell_candidate(graph, candidates[candidate_id]).sequence
        digest = _sha256_sequence(sequence)
        if digest in sequence_by_digest and sequence_by_digest[digest] != sequence:
            raise RuntimeError("SHA256 collision in candidate sequence deduplication")
        sequence_by_digest[digest] = sequence
        representative_by_digest.setdefault(digest, candidate_id)
        digest_by_candidate[candidate_id] = digest

    class_id_by_digest = {
        digest: f"BGCPHASER_SEQ{number:06d}"
        for number, digest in enumerate(sorted(representative_by_digest), start=1)
    }
    class_by_candidate = {
        candidate_id: class_id_by_digest[digest]
        for candidate_id, digest in digest_by_candidate.items()
    }
    representative_by_class = {
        class_id_by_digest[digest]: candidate_id
        for digest, candidate_id in representative_by_digest.items()
    }
    digest_by_class = {
        class_id_by_digest[digest]: digest
        for digest in representative_by_digest
    }
    length_by_class = {
        class_id_by_digest[digest]: len(sequence_by_digest[digest])
        for digest in representative_by_digest
    }
    return (
        class_by_candidate,
        representative_by_class,
        digest_by_class,
        length_by_class,
    )


def _write_unique_candidate_fasta(
    path: Path,
    *,
    graph: nx.DiGraph,
    candidates: Mapping[str, Sequence[OrientedNode]],
    representative_by_class: Mapping[str, str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for class_id in sorted(representative_by_class):
            representative = representative_by_class[class_id]
            sequence = spell_candidate(
                graph,
                candidates[representative],
            ).sequence
            handle.write(f">{class_id}\n{sequence}\n")


def _write_transition_backmap(path: Path, classes: TransitionClasses) -> None:
    rows = []
    for occurrence in classes.occurrences:
        rows.append(
            {
                "candidate_id": occurrence.candidate_id,
                "transition_index": occurrence.transition_index,
                "left_state": occurrence.left_state[0] + occurrence.left_state[1],
                "right_state": occurrence.right_state[0] + occurrence.right_state[1],
                "R_class_id": occurrence.r_class_id or "NA",
                "L_class_id": occurrence.l_class_id or "NA",
            }
        )
    _write_tsv(
        path,
        [
            "candidate_id",
            "transition_index",
            "left_state",
            "right_state",
            "R_class_id",
            "L_class_id",
        ],
        rows,
    )


def analyse_candidates(
    *,
    gfa_path: str | Path,
    candidate_path: str | Path,
    read1_path: str | Path,
    read2_path: str | Path,
    output_dir: str | Path,
    insert_min: int,
    insert_max: int,
    antismash_database_root: str | Path,
    unpaired_read_paths: Sequence[str | Path] = (),
    assembler_paths: Sequence[str | Path] = (),
    depth_tag: str | None = None,
    minimap2_executable: str | Path = "minimap2",
    antismash_executable: str | Path = "antismash",
    minimap2_threads: int = 1,
    antismash_cpus: int = 1,
) -> AnalysisRun:
    """Execute the frozen BGCPhaser intrinsic scoring pipeline.

    Reference or truth sequences are intentionally not accepted as inputs.
    """
    if insert_min < 0:
        raise ValueError("insert_min must be >= 0")
    if insert_max < insert_min:
        raise ValueError("insert_max must be >= insert_min")
    if minimap2_threads < 1:
        raise ValueError("minimap2_threads must be >= 1")
    if antismash_cpus < 1:
        raise ValueError("antismash_cpus must be >= 1")

    gfa = _require_file(gfa_path, label="GFA")
    paired_reads = [
        _require_file(read1_path, label="read1"),
        _require_file(read2_path, label="read2"),
    ]
    unpaired_reads = [
        _require_file(path, label=f"unpaired_read_{index}")
        for index, path in enumerate(unpaired_read_paths, start=1)
    ]
    database_root = Path(antismash_database_root)
    if not database_root.is_dir():
        raise FileNotFoundError(
            f"antiSMASH database directory absent: {database_root}"
        )

    graph = read_gfa(gfa)
    candidates = _read_candidates(candidate_path)
    for candidate_id in sorted(candidates):
        try:
            spell_candidate(graph, candidates[candidate_id])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"{candidate_id}: {exc}") from exc

    depth_values = _depth_by_segment(graph, depth_tag)
    assembler_adjacency = _assembler_adjacency(assembler_paths)
    transition_classes = construct_transition_classes(
        graph,
        candidates,
        target_flank=insert_max,
        minimum_side=R_CONTEXT_BP,
        class_prefix="BGCPHASER",
    )
    (
        sequence_class_by_candidate,
        representative_by_class,
        digest_by_class,
        sequence_length_by_class,
    ) = _candidate_sequence_classes(graph, candidates)

    final_dir = Path(output_dir)
    staging_dir = Path(str(final_dir) + ".staging")
    if final_dir.exists():
        raise FileExistsError(f"Output directory already exists: {final_dir}")
    if staging_dir.exists():
        raise FileExistsError(
            f"Staging directory already exists from an earlier attempt: {staging_dir}"
        )
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir.mkdir()
    evidence_dir = staging_dir / "evidence"
    chemistry_dir = staging_dir / "chemistry"
    candidate_sequence_dir = staging_dir / "candidate_sequences"
    evidence_dir.mkdir()
    chemistry_dir.mkdir()
    candidate_sequence_dir.mkdir()

    r_fasta = evidence_dir / "r_classes.fasta"
    l_fasta = evidence_dir / "l_classes.fasta"
    write_fasta(transition_classes.r_sequences, r_fasta)
    write_fasta(transition_classes.l_sequences, l_fasta)
    _write_transition_backmap(
        evidence_dir / "transition_class_backmap.tsv",
        transition_classes,
    )

    r_support: Counter[str] = Counter()
    if transition_classes.r_sequences:
        r_target_lengths = {
            class_id: len(sequence)
            for class_id, sequence in transition_classes.r_sequences.items()
        }
        r_sources = [*paired_reads, *unpaired_reads]
        for source_index, read_path in enumerate(r_sources, start=1):
            r_sam = evidence_dir / f"r_source_{source_index:03d}.sam"
            run_minimap2(
                executable=minimap2_executable,
                target_fasta=r_fasta,
                reads=[read_path],
                sam_path=r_sam,
                log_path=(
                    evidence_dir
                    / f"r_source_{source_index:03d}_minimap2.log"
                ),
                threads=minimap2_threads,
            )
            r_support.update(
                _r_support_from_sam(r_sam, r_target_lengths)
            )

    l_support: Counter[str] = Counter()
    if transition_classes.l_sequences:
        l_sam = evidence_dir / "l_alignments.sam"
        run_minimap2(
            executable=minimap2_executable,
            target_fasta=l_fasta,
            reads=paired_reads,
            sam_path=l_sam,
            log_path=evidence_dir / "l_minimap2.log",
            threads=minimap2_threads,
        )
        l_support = _l_support_from_sam(
            l_sam,
            junction_by_class=transition_classes.l_junctions,
            insert_min=insert_min,
            insert_max=insert_max,
        )

    _write_tsv(
        evidence_dir / "r_class_support.tsv",
        ["class_id", "sequence_length_bp", "supporting_read_end_count"],
        [
            {
                "class_id": class_id,
                "sequence_length_bp": len(
                    transition_classes.r_sequences[class_id]
                ),
                "supporting_read_end_count": int(
                    r_support.get(class_id, 0)
                ),
            }
            for class_id in sorted(transition_classes.r_sequences)
        ],
    )
    _write_tsv(
        evidence_dir / "l_class_support.tsv",
        [
            "class_id",
            "sequence_length_bp",
            "junction_0based",
            "supporting_fragment_count",
        ],
        [
            {
                "class_id": class_id,
                "sequence_length_bp": len(
                    transition_classes.l_sequences[class_id]
                ),
                "junction_0based": transition_classes.l_junctions[class_id],
                "supporting_fragment_count": int(
                    l_support.get(class_id, 0)
                ),
            }
            for class_id in sorted(transition_classes.l_sequences)
        ],
    )

    unique_fasta = candidate_sequence_dir / "unique_candidate_sequences.fasta"
    _write_unique_candidate_fasta(
        unique_fasta,
        graph=graph,
        candidates=candidates,
        representative_by_class=representative_by_class,
    )
    _write_tsv(
        candidate_sequence_dir / "candidate_sequence_backmap.tsv",
        [
            "candidate_id",
            "sequence_class_id",
            "sequence_sha256",
            "sequence_length_bp",
        ],
        [
            {
                "candidate_id": candidate_id,
                "sequence_class_id": sequence_class_by_candidate[candidate_id],
                "sequence_sha256": digest_by_class[
                    sequence_class_by_candidate[candidate_id]
                ],
                "sequence_length_bp": sequence_length_by_class[
                    sequence_class_by_candidate[candidate_id]
                ],
            }
            for candidate_id in sorted(candidates)
        ],
    )

    chemistry_by_class: dict[str, ChemistryFeatures] = {}
    for sequence_class_id in sorted(representative_by_class):
        representative = representative_by_class[sequence_class_id]
        sequence = spell_candidate(
            graph,
            candidates[representative],
        ).sequence
        input_fasta = (
            candidate_sequence_dir
            / "antismash_inputs"
            / f"{sequence_class_id}.fasta"
        )
        write_fasta({sequence_class_id: sequence}, input_fasta)
        json_path = run_antismash(
            executable=antismash_executable,
            input_fasta=input_fasta,
            output_dir=chemistry_dir / sequence_class_id,
            output_basename=sequence_class_id,
            database_root=database_root,
            cpus=antismash_cpus,
        )
        chemistry_by_class[sequence_class_id] = (
            chemistry_from_antismash_json(json_path)
        )

    _write_tsv(
        chemistry_dir / "chemistry_classes.tsv",
        [
            "sequence_class_id",
            "sequence_sha256",
            "C",
            "C_status",
            "C_supported_module_count",
            "C_generic_unresolved_module_count",
            "C_prediction_missing_module_count",
            "C_assessable_explicit_loader_module_count",
            "M",
            "M_status",
            "M_complete_module_count",
            "M_incomplete_module_count",
            "M_total_module_count",
        ],
        [
            {
                "sequence_class_id": class_id,
                "sequence_sha256": digest_by_class[class_id],
                "C": format_q12(chemistry_by_class[class_id].C),
                "C_status": chemistry_by_class[class_id].C_status,
                "C_supported_module_count": chemistry_by_class[
                    class_id
                ].C_supported_module_count,
                "C_generic_unresolved_module_count": chemistry_by_class[
                    class_id
                ].C_generic_unresolved_module_count,
                "C_prediction_missing_module_count": chemistry_by_class[
                    class_id
                ].C_prediction_missing_module_count,
                "C_assessable_explicit_loader_module_count": chemistry_by_class[
                    class_id
                ].C_assessable_explicit_loader_module_count,
                "M": format_q12(chemistry_by_class[class_id].M),
                "M_status": chemistry_by_class[class_id].M_status,
                "M_complete_module_count": chemistry_by_class[
                    class_id
                ].M_complete_module_count,
                "M_incomplete_module_count": chemistry_by_class[
                    class_id
                ].M_incomplete_module_count,
                "M_total_module_count": chemistry_by_class[
                    class_id
                ].M_total_module_count,
            }
            for class_id in sorted(chemistry_by_class)
        ],
    )

    sequence_detail: dict[
        str,
        tuple[
            object,
            object,
            CoverageConsistency,
            AssemblerDiscordance,
            SequenceComposite,
        ],
    ] = {}
    chemistry_detail: dict[str, ChemistryFeatures] = {}
    scored = []

    for candidate_id in sorted(candidates):
        r_value = candidate_transition_support(
            transition_classes.occurrences,
            candidate_id=candidate_id,
            evidence="R",
            class_support=r_support,
        )
        l_value = candidate_transition_support(
            transition_classes.occurrences,
            candidate_id=candidate_id,
            evidence="L",
            class_support=l_support,
        )
        coverage = coverage_consistency(
            candidates[candidate_id],
            depth_values,
        )
        discordance = assembler_discordance(
            candidates[candidate_id],
            assembler_adjacency,
        )
        sequence_score = sequence_composite(
            r_value,
            l_value,
            coverage.V,
            discordance.H,
        )
        chemistry = chemistry_by_class[
            sequence_class_by_candidate[candidate_id]
        ]
        sequence_detail[candidate_id] = (
            r_value,
            l_value,
            coverage,
            discordance,
            sequence_score,
        )
        chemistry_detail[candidate_id] = chemistry
        scored.append(
            score_candidate(
                candidate_id,
                sequence_score.score,
                chemistry.C,
                chemistry.M,
            )
        )

    ranked = rank_candidates(scored)
    ranked_path = staging_dir / "ranked_candidates.tsv"
    ranked_rows = []
    for item in ranked:
        r_value, l_value, coverage, discordance, sequence_score = (
            sequence_detail[item.candidate_id]
        )
        chemistry = chemistry_detail[item.candidate_id]
        sequence_class_id = sequence_class_by_candidate[item.candidate_id]
        ranked_rows.append(
            {
                "rank": item.rank if item.rank is not None else "NA",
                "candidate_id": item.candidate_id,
                "oriented_walk": format_oriented_path(
                    candidates[item.candidate_id]
                ),
                "sequence_class_id": sequence_class_id,
                "sequence_sha256": digest_by_class[sequence_class_id],
                "R": format_q12(r_value),
                "L": format_q12(l_value),
                "V": format_q12(coverage.V),
                "H": format_q12(discordance.H),
                "one_minus_H": format_q12(sequence_score.one_minus_H),
                "sequence_defined_components": ";".join(
                    sequence_score.defined_components
                ),
                "sequence_score": format_q12(item.sequence_score),
                "chemistry_c": format_q12(item.chemistry_c),
                "chemistry_c_status": chemistry.C_status,
                "chemistry_m": format_q12(item.chemistry_m),
                "chemistry_m_status": chemistry.M_status,
                "chemistry_score": format_q12(item.chemistry_score),
                "combined_score": format_q12(item.combined_score),
                "status": item.status,
            }
        )

    _write_tsv(
        ranked_path,
        [
            "rank",
            "candidate_id",
            "oriented_walk",
            "sequence_class_id",
            "sequence_sha256",
            "R",
            "L",
            "V",
            "H",
            "one_minus_H",
            "sequence_defined_components",
            "sequence_score",
            "chemistry_c",
            "chemistry_c_status",
            "chemistry_m",
            "chemistry_m_status",
            "chemistry_score",
            "combined_score",
            "status",
        ],
        ranked_rows,
    )

    _write_tsv(
        staging_dir / "run_parameters.tsv",
        ["parameter", "value"],
        [
            {"parameter": "R_context_bp", "value": R_CONTEXT_BP},
            {"parameter": "insert_min_bp", "value": insert_min},
            {"parameter": "insert_max_bp", "value": insert_max},
            {"parameter": "L_target_flank_bp", "value": insert_max},
            {
                "parameter": "depth_tag",
                "value": depth_tag if depth_tag is not None else "NA",
            },
            {
                "parameter": "assembler_path_file_count",
                "value": len(assembler_paths),
            },
            {
                "parameter": "unpaired_R_read_file_count",
                "value": len(unpaired_reads),
            },
            {
                "parameter": "R_read_source_count",
                "value": len(paired_reads) + len(unpaired_reads),
            },
            {"parameter": "reference_input_used", "value": "NO"},
            {
                "parameter": "reference_comparison_used_in_scoring",
                "value": "NO",
            },
        ],
    )

    ranked_count = sum(item.rank is not None for item in ranked)
    os.replace(staging_dir, final_dir)
    return AnalysisRun(
        output_dir=final_dir,
        ranked_path=final_dir / "ranked_candidates.tsv",
        candidate_count=len(candidates),
        ranked_count=ranked_count,
    )
