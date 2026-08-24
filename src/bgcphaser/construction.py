from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import hashlib
import re
from pathlib import Path
from typing import Mapping, Sequence

import networkx as nx

from .gfa import (
    OrientedNode,
    OrientedPath,
    reverse_complement,
)
from .sequence import support_fraction


_OVERLAP_RE = re.compile(
    r"([0-9]+)M"
)


@dataclass(frozen=True)
class SpelledCandidate:
    sequence: str
    junctions: tuple[int, ...]
    overlaps: tuple[int, ...]


@dataclass(frozen=True)
class TransitionOccurrence:
    candidate_id: str
    transition_index: int
    left_state: OrientedNode
    right_state: OrientedNode
    r_class_id: str | None
    l_class_id: str | None


@dataclass(frozen=True)
class TransitionClasses:
    r_sequences: dict[str, str]
    l_sequences: dict[str, str]
    l_junctions: dict[str, int]
    occurrences: tuple[
        TransitionOccurrence,
        ...,
    ]


def parse_oriented_path(
    text: str,
) -> OrientedPath:
    """
    Parse the frozen BGCPhaser oriented-path notation:

        segment+->segment-->segment+
    """
    path: OrientedPath = []

    for token in text.split(
        "->"
    ):
        if len(token) < 2:
            raise ValueError(
                f"Malformed oriented state: "
                f"{token!r}"
            )

        orientation = token[-1]
        segment_id = token[:-1]

        if orientation not in {
            "+",
            "-",
        }:
            raise ValueError(
                f"Invalid orientation: "
                f"{token!r}"
            )

        if not segment_id:
            raise ValueError(
                "Missing segment ID"
            )

        path.append(
            (
                segment_id,
                orientation,
            )
        )

    if not path:
        raise ValueError(
            "Empty oriented path"
        )

    return path


def _overlap_bp(
    graph: nx.DiGraph,
    left: OrientedNode,
    right: OrientedNode,
) -> int:
    if not graph.has_edge(
        left,
        right,
    ):
        raise ValueError(
            "Candidate contains missing "
            f"GFA edge {left}->{right}"
        )

    attributes = graph.edges[
        left,
        right,
    ]

    recorded_overlaps = (
        attributes.get(
            "overlaps"
        )
    )

    if recorded_overlaps is not None:
        distinct_overlaps = set(
            recorded_overlaps
        )

        if len(
            distinct_overlaps
        ) != 1:
            raise ValueError(
                "Conflicting overlap CIGARs "
                f"for GFA edge "
                f"{left}->{right}: "
                f"{sorted(distinct_overlaps)}"
            )

    overlap = attributes.get(
        "overlap"
    )

    if not isinstance(
        overlap,
        str,
    ):
        raise ValueError(
            "GFA transition lacks a "
            "defined overlap CIGAR"
        )

    match = _OVERLAP_RE.fullmatch(
        overlap
    )

    if match is None:
        raise ValueError(
            "Unsupported GFA overlap "
            f"CIGAR: {overlap}"
        )

    return int(
        match.group(1)
    )


def spell_candidate(
    graph: nx.DiGraph,
    path: Sequence[OrientedNode],
) -> SpelledCandidate:
    """
    Spell one candidate using the validated coordinate
    convention.

    For every transition, the junction coordinate is the
    current spelled length before adding the non-overlapping
    suffix of the right segment.
    """
    if not path:
        raise ValueError(
            "Cannot spell an empty path"
        )

    for state in path:
        if state not in graph:
            raise KeyError(
                f"Oriented state absent "
                f"from GFA: {state}"
            )

    first_sequence = graph.nodes[
        path[0]
    ].get(
        "sequence"
    )

    if not isinstance(
        first_sequence,
        str,
    ):
        raise ValueError(
            "GFA segment sequence absent"
        )

    sequence = first_sequence
    junctions: list[int] = []
    overlaps: list[int] = []

    for left, right in zip(
        path,
        path[1:],
    ):
        overlap = _overlap_bp(
            graph,
            left,
            right,
        )

        right_sequence = graph.nodes[
            right
        ].get(
            "sequence"
        )

        if not isinstance(
            right_sequence,
            str,
        ):
            raise ValueError(
                "GFA segment sequence absent"
            )

        if overlap > len(
            right_sequence
        ):
            raise ValueError(
                "GFA overlap exceeds "
                "right-segment length"
            )

        junctions.append(
            len(sequence)
        )

        overlaps.append(
            overlap
        )

        sequence += (
            right_sequence[
                overlap:
            ]
        )

    return SpelledCandidate(
        sequence=sequence,
        junctions=tuple(
            junctions
        ),
        overlaps=tuple(
            overlaps
        ),
    )


def _sha256_sequence(
    sequence: str,
) -> str:
    return hashlib.sha256(
        sequence.encode(
            "ascii"
        )
    ).hexdigest()


def canonical_sequence(
    sequence: str,
) -> str:
    """
    Strand-canonical nucleotide representation.
    """
    forward = sequence.upper()
    reverse = reverse_complement(
        forward
    ).upper()

    return min(
        forward,
        reverse,
    )


def canonical_sequence_with_junction(
    sequence: str,
    junction: int,
) -> tuple[str, int]:
    """
    Strand-canonicalise a transition target while
    preserving the junction boundary coordinate.
    """
    if not (
        0 <= junction <= len(sequence)
    ):
        raise ValueError(
            "Junction outside sequence"
        )

    forward = sequence.upper()
    reverse = reverse_complement(
        forward
    ).upper()

    if forward <= reverse:
        return (
            forward,
            junction,
        )

    return (
        reverse,
        len(sequence) - junction,
    )


def construct_transition_classes(
    graph: nx.DiGraph,
    candidates: Mapping[
        str,
        Sequence[OrientedNode],
    ],
    *,
    target_flank: int,
    minimum_side: int = 30,
    class_prefix: str = "BGCPHASER",
) -> TransitionClasses:
    """
    Construct the validated R/L primitive evidence classes.

    R:
      junction witness =
      [junction - overlap - minimum_side,
       junction + minimum_side]

    L:
      candidate transition target =
      up to target_flank bp on each side of the
      junction, then strand-canonicalised together
      with its junction coordinate.

    Classes are deduplicated by SHA256 of the
    canonical nucleotide sequence and numbered in
    sorted-digest order.
    """
    if target_flank < minimum_side:
        raise ValueError(
            "target_flank must be >= "
            "minimum_side"
        )

    if minimum_side < 1:
        raise ValueError(
            "minimum_side must be >= 1"
        )

    r_by_digest: dict[
        str,
        str,
    ] = {}

    l_by_digest: dict[
        str,
        str,
    ] = {}

    l_junction_by_digest: dict[
        str,
        int,
    ] = {}

    raw_occurrences: list[
        tuple[
            str,
            int,
            OrientedNode,
            OrientedNode,
            str | None,
            str | None,
        ]
    ] = []

    for candidate_id in sorted(
        candidates
    ):
        path = list(
            candidates[
                candidate_id
            ]
        )

        spelling = spell_candidate(
            graph,
            path,
        )

        transitions = list(
            zip(
                path,
                path[1:],
            )
        )

        if not (
            len(transitions)
            == len(
                spelling.junctions
            )
            == len(
                spelling.overlaps
            )
        ):
            raise RuntimeError(
                "Transition/spelling "
                "cardinality drift"
            )

        for (
            transition_index,
            (
                left,
                right,
            ),
            junction,
            overlap,
        ) in zip(
            range(
                1,
                len(transitions) + 1,
            ),
            transitions,
            spelling.junctions,
            spelling.overlaps,
        ):
            target_start = max(
                0,
                junction - target_flank,
            )

            target_end = min(
                len(
                    spelling.sequence
                ),
                junction + target_flank,
            )

            target_junction = (
                junction
                - target_start
            )

            left_available = (
                target_junction
            )

            right_available = (
                target_end
                - junction
            )

            assessable = (
                left_available
                >= minimum_side
                and right_available
                >= minimum_side
            )

            if not assessable:
                raw_occurrences.append(
                    (
                        candidate_id,
                        transition_index,
                        left,
                        right,
                        None,
                        None,
                    )
                )

                continue

            target_sequence = (
                spelling.sequence[
                    target_start:
                    target_end
                ]
            )

            (
                canonical_l,
                canonical_l_junction,
            ) = (
                canonical_sequence_with_junction(
                    target_sequence,
                    target_junction,
                )
            )

            l_digest = _sha256_sequence(
                canonical_l
            )

            previous_l = l_by_digest.get(
                l_digest
            )

            if (
                previous_l is not None
                and previous_l
                != canonical_l
            ):
                raise RuntimeError(
                    "SHA256 collision in "
                    "L class construction"
                )

            previous_junction = (
                l_junction_by_digest.get(
                    l_digest
                )
            )

            if (
                previous_junction
                is not None
                and previous_junction
                != canonical_l_junction
            ):
                raise ValueError(
                    "Identical L class sequence "
                    "has inconsistent junction "
                    "coordinate"
                )

            l_by_digest[
                l_digest
            ] = canonical_l

            l_junction_by_digest[
                l_digest
            ] = canonical_l_junction

            r_start = (
                junction
                - overlap
                - minimum_side
            )

            r_end = (
                junction
                + minimum_side
            )

            if (
                r_start < 0
                or r_end
                > len(
                    spelling.sequence
                )
            ):
                raise ValueError(
                    f"{candidate_id}: "
                    "R witness not constructible "
                    "for an otherwise assessable "
                    "transition"
                )

            r_sequence = (
                spelling.sequence[
                    r_start:
                    r_end
                ]
            )

            expected_r_length = (
                overlap
                + 2 * minimum_side
            )

            if (
                len(r_sequence)
                != expected_r_length
            ):
                raise RuntimeError(
                    "R witness length drift"
                )

            canonical_r = (
                canonical_sequence(
                    r_sequence
                )
            )

            r_digest = _sha256_sequence(
                canonical_r
            )

            previous_r = r_by_digest.get(
                r_digest
            )

            if (
                previous_r is not None
                and previous_r
                != canonical_r
            ):
                raise RuntimeError(
                    "SHA256 collision in "
                    "R class construction"
                )

            r_by_digest[
                r_digest
            ] = canonical_r

            raw_occurrences.append(
                (
                    candidate_id,
                    transition_index,
                    left,
                    right,
                    r_digest,
                    l_digest,
                )
            )

    r_id_by_digest = {
        digest: (
            f"{class_prefix}_RCL"
            f"{number:04d}"
        )
        for number, digest
        in enumerate(
            sorted(
                r_by_digest
            ),
            start=1,
        )
    }

    l_id_by_digest = {
        digest: (
            f"{class_prefix}_LCL"
            f"{number:04d}"
        )
        for number, digest
        in enumerate(
            sorted(
                l_by_digest
            ),
            start=1,
        )
    }

    r_sequences = {
        r_id_by_digest[
            digest
        ]:
            r_by_digest[
                digest
            ]
        for digest in sorted(
            r_by_digest
        )
    }

    l_sequences = {
        l_id_by_digest[
            digest
        ]:
            l_by_digest[
                digest
            ]
        for digest in sorted(
            l_by_digest
        )
    }

    l_junctions = {
        l_id_by_digest[
            digest
        ]:
            l_junction_by_digest[
                digest
            ]
        for digest in sorted(
            l_by_digest
        )
    }

    occurrences = tuple(
        TransitionOccurrence(
            candidate_id=(
                candidate_id
            ),
            transition_index=(
                transition_index
            ),
            left_state=left,
            right_state=right,
            r_class_id=(
                r_id_by_digest[
                    r_digest
                ]
                if r_digest is not None
                else None
            ),
            l_class_id=(
                l_id_by_digest[
                    l_digest
                ]
                if l_digest is not None
                else None
            ),
        )
        for (
            candidate_id,
            transition_index,
            left,
            right,
            r_digest,
            l_digest,
        ) in raw_occurrences
    )

    return TransitionClasses(
        r_sequences=r_sequences,
        l_sequences=l_sequences,
        l_junctions=l_junctions,
        occurrences=occurrences,
    )


def write_fasta(
    records: Mapping[str, str],
    path: str | Path,
) -> Path:
    """
    Write deterministic single-line FASTA records.
    """
    output = Path(path)

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output.open(
        "w",
        encoding="ascii",
        newline="\n",
    ) as handle:
        for record_id in sorted(
            records
        ):
            sequence = records[
                record_id
            ]

            handle.write(
                f">{record_id}\n"
                f"{sequence}\n"
            )

    return output


def candidate_transition_support(
    occurrences: Sequence[
        TransitionOccurrence
    ],
    *,
    candidate_id: str,
    evidence: str,
    class_support: Mapping[
        str,
        int,
    ],
) -> Decimal | None:
    """
    Propagate primitive class support to one candidate.

    Support-count magnitude is not used. Each assessable
    transition occurrence is binary:
      class support count > 0 -> supported
      class support count = 0 -> unsupported
    """
    if evidence not in {
        "R",
        "L",
    }:
        raise ValueError(
            "evidence must be 'R' or 'L'"
        )

    class_ids: list[str] = []

    for occurrence in occurrences:
        if (
            occurrence.candidate_id
            != candidate_id
        ):
            continue

        class_id = (
            occurrence.r_class_id
            if evidence == "R"
            else occurrence.l_class_id
        )

        if class_id is not None:
            class_ids.append(
                class_id
            )

    if not class_ids:
        return None

    supported = sum(
        class_support.get(
            class_id,
            0,
        ) > 0
        for class_id in class_ids
    )

    return support_fraction(
        supported,
        len(class_ids),
    )
