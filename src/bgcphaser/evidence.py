from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import re
from typing import Iterable, Mapping, Sequence


_CIGAR_RE = re.compile(
    r"([0-9]+)([MIDNSHP=X])"
)


@dataclass(frozen=True)
class Alignment:
    qname: str
    target: str
    cigar: str
    alignment_score: int | None
    start0: int | None
    end0: int | None
    unmapped: bool
    reverse: bool
    secondary: bool
    supplementary: bool
    first: bool
    second: bool


@dataclass(frozen=True)
class SupportDecision:
    decision: str
    supported_class_id: str | None
    best_score: int
    second_best_score: int | None
    score_margin: int | None
    competing_class_count: int
    tied_best_class_ids: tuple[str, ...]


@dataclass(frozen=True)
class FragmentGeometry:
    score: int
    fragment_length: int
    outer_start0: int
    outer_end0: int

    mate1_start0: int
    mate1_end0: int
    mate1_reverse: bool
    mate1_score: int

    mate2_start0: int
    mate2_end0: int
    mate2_reverse: bool
    mate2_score: int


def parse_cigar(
    cigar: str,
) -> list[tuple[int, str]]:
    if cigar == "*":
        return []

    pieces = [
        (
            int(length),
            operation,
        )
        for length, operation
        in _CIGAR_RE.findall(cigar)
    ]

    rebuilt = "".join(
        f"{length}{operation}"
        for length, operation in pieces
    )

    if rebuilt != cigar:
        raise ValueError(
            f"Cannot parse CIGAR: {cigar}"
        )

    return pieces


def reference_end(
    start0: int,
    cigar: str,
) -> int:
    position = start0

    for length, operation in parse_cigar(
        cigar
    ):
        if operation in {
            "M",
            "D",
            "N",
            "=",
            "X",
        }:
            position += length

    return position


def parse_sam_record(
    raw: str,
) -> Alignment:
    fields = raw.rstrip(
        "\n"
    ).split(
        "\t"
    )

    if len(fields) < 11:
        raise ValueError(
            "Malformed SAM record"
        )

    flag = int(fields[1])

    tags: dict[
        str,
        tuple[str, str],
    ] = {}

    for token in fields[11:]:
        pieces = token.split(
            ":",
            2,
        )

        if len(pieces) == 3:
            tags[pieces[0]] = (
                pieces[1],
                pieces[2],
            )

    unmapped = bool(
        flag & 0x4
    )

    supplementary = bool(
        flag & 0x800
    )

    if unmapped:
        start0 = None
        end0 = None

    else:
        start0 = (
            int(fields[3]) - 1
        )

        end0 = reference_end(
            start0,
            fields[5],
        )

    score: int | None = None

    if "AS" in tags:
        tag_type, raw_score = tags["AS"]

        if tag_type != "i":
            raise ValueError(
                "Unexpected SAM AS tag type"
            )

        score = int(raw_score)

    return Alignment(
        qname=fields[0],
        target=fields[2],
        cigar=fields[5],
        alignment_score=score,
        start0=start0,
        end0=end0,
        unmapped=unmapped,
        reverse=bool(flag & 0x10),
        secondary=bool(flag & 0x100),
        supplementary=supplementary,
        first=bool(flag & 0x40),
        second=bool(flag & 0x80),
    )


def covers_full_target(
    alignment: Alignment,
    target_length: int,
) -> bool:
    """
    Validated R witness-coverage rule.

    The alignment must begin at reference coordinate 0
    and aligned M/= /X operations must collectively cover
    every target position.

    Deletions and reference skips invalidate the witness.
    Insertions/clipping do not advance the target position.
    """
    if (
        alignment.unmapped
        or alignment.start0 != 0
    ):
        return False

    position = 0

    for length, operation in parse_cigar(
        alignment.cigar
    ):
        if operation in {
            "M",
            "=",
            "X",
        }:
            position += length

        elif operation in {
            "D",
            "N",
        }:
            return False

        elif operation in {
            "I",
            "S",
            "H",
            "P",
        }:
            continue

        else:
            raise ValueError(
                f"Unexpected CIGAR op: {operation}"
            )

    return position == target_length


def strict_best_class(
    best_by_class: Mapping[str, int],
) -> SupportDecision | None:
    """
    Frozen strict-AS competition.

    One uniquely highest class receives support.
    Exact best-score ties are ambiguous.
    No minimum score margin is imposed.
    """
    if not best_by_class:
        return None

    best_score = max(
        best_by_class.values()
    )

    winners = tuple(
        sorted(
            class_id
            for class_id, score
            in best_by_class.items()
            if score == best_score
        )
    )

    ordered_scores = sorted(
        best_by_class.values(),
        reverse=True,
    )

    second_score = (
        ordered_scores[1]
        if len(ordered_scores) > 1
        else None
    )

    if len(winners) == 1:
        margin = (
            best_score - second_score
            if second_score is not None
            else None
        )

        return SupportDecision(
            decision=(
                "STRICT_UNIQUE_BEST_SUPPORT"
            ),
            supported_class_id=winners[0],
            best_score=best_score,
            second_best_score=second_score,
            score_margin=margin,
            competing_class_count=len(
                best_by_class
            ),
            tied_best_class_ids=(),
        )

    return SupportDecision(
        decision="TIED_BEST_AMBIGUOUS",
        supported_class_id=None,
        best_score=best_score,
        second_best_score=best_score,
        score_margin=0,
        competing_class_count=len(
            best_by_class
        ),
        tied_best_class_ids=winners,
    )


def r_support_decision(
    alignments: Iterable[Alignment],
    target_lengths: Mapping[str, int],
) -> SupportDecision | None:
    """
    Reduce one read end to the validated R support decision.
    """
    best_by_class: dict[
        str,
        int,
    ] = {}

    for alignment in alignments:
        if (
            alignment.unmapped
            or alignment.supplementary
        ):
            continue

        if alignment.target not in target_lengths:
            raise ValueError(
                "Unexpected R class target: "
                f"{alignment.target}"
            )

        if alignment.alignment_score is None:
            raise ValueError(
                "Mapped non-supplementary "
                "R alignment lacks AS"
            )

        if not covers_full_target(
            alignment,
            target_lengths[
                alignment.target
            ],
        ):
            continue

        previous = best_by_class.get(
            alignment.target
        )

        if (
            previous is None
            or alignment.alignment_score
            > previous
        ):
            best_by_class[
                alignment.target
            ] = alignment.alignment_score

    return strict_best_class(
        best_by_class
    )


def inward_fr(
    first: Alignment,
    second: Alignment,
) -> bool:
    """
    Frozen inward-facing FR predicate.
    """
    if (
        first.start0 is None
        or second.start0 is None
    ):
        return False

    if (
        first.reverse
        == second.reverse
    ):
        return False

    if not first.reverse:
        forward = first
        reverse = second
    else:
        forward = second
        reverse = first

    return (
        forward.start0
        < reverse.start0
    )


def best_geometry_for_class(
    mate1: Sequence[Alignment],
    mate2: Sequence[Alignment],
    *,
    junction: int,
    insert_min: int,
    insert_max: int,
) -> FragmentGeometry | None:
    """
    Frozen L geometry and within-class maximisation.

    Valid pair:
      - inward-facing FR
      - outer span crosses junction
      - outer span within configured insert interval

    Class score = maximum AS1 + AS2.
    """
    if insert_min < 0:
        raise ValueError(
            "insert_min must be >= 0"
        )

    if insert_max < insert_min:
        raise ValueError(
            "insert_max must be >= insert_min"
        )

    best: FragmentGeometry | None = None

    for first in mate1:
        for second in mate2:
            if (
                first.unmapped
                or second.unmapped
                or first.supplementary
                or second.supplementary
            ):
                continue

            if (
                first.start0 is None
                or first.end0 is None
                or second.start0 is None
                or second.end0 is None
            ):
                continue

            if (
                first.alignment_score is None
                or second.alignment_score is None
            ):
                raise ValueError(
                    "Mapped non-supplementary "
                    "L alignment lacks AS"
                )

            if not inward_fr(
                first,
                second,
            ):
                continue

            outer_start = min(
                first.start0,
                second.start0,
            )

            outer_end = max(
                first.end0,
                second.end0,
            )

            fragment_length = (
                outer_end
                - outer_start
            )

            if not (
                insert_min
                <= fragment_length
                <= insert_max
            ):
                continue

            if not (
                outer_start
                < junction
                < outer_end
            ):
                continue

            score = (
                first.alignment_score
                + second.alignment_score
            )

            candidate = FragmentGeometry(
                score=score,
                fragment_length=fragment_length,
                outer_start0=outer_start,
                outer_end0=outer_end,
                mate1_start0=first.start0,
                mate1_end0=first.end0,
                mate1_reverse=first.reverse,
                mate1_score=(
                    first.alignment_score
                ),
                mate2_start0=second.start0,
                mate2_end0=second.end0,
                mate2_reverse=second.reverse,
                mate2_score=(
                    second.alignment_score
                ),
            )

            if (
                best is None
                or candidate.score
                > best.score
            ):
                best = candidate

            elif (
                candidate.score
                == best.score
            ):
                old_key = (
                    best.fragment_length,
                    best.outer_start0,
                    best.outer_end0,
                    best.mate1_start0,
                    best.mate2_start0,
                )

                new_key = (
                    candidate.fragment_length,
                    candidate.outer_start0,
                    candidate.outer_end0,
                    candidate.mate1_start0,
                    candidate.mate2_start0,
                )

                if new_key < old_key:
                    best = candidate

    return best


def l_support_decision(
    alignments: Iterable[Alignment],
    *,
    junction_by_class: Mapping[str, int],
    insert_min: int,
    insert_max: int,
) -> tuple[
    SupportDecision | None,
    FragmentGeometry | None,
]:
    """
    Reduce one paired fragment to the validated L decision.
    """
    usable = [
        alignment
        for alignment in alignments
        if (
            not alignment.unmapped
            and not alignment.supplementary
        )
    ]

    for alignment in usable:
        if (
            alignment.target
            not in junction_by_class
        ):
            raise ValueError(
                "Unexpected L class target: "
                f"{alignment.target}"
            )

        if alignment.alignment_score is None:
            raise ValueError(
                "Mapped non-supplementary "
                "L alignment lacks AS"
            )

    by_target_1: dict[
        str,
        list[Alignment],
    ] = defaultdict(list)

    by_target_2: dict[
        str,
        list[Alignment],
    ] = defaultdict(list)

    for alignment in usable:
        if alignment.first:
            by_target_1[
                alignment.target
            ].append(alignment)

        if alignment.second:
            by_target_2[
                alignment.target
            ].append(alignment)

    class_best: dict[
        str,
        FragmentGeometry,
    ] = {}

    for class_id in (
        set(by_target_1)
        & set(by_target_2)
    ):
        geometry = best_geometry_for_class(
            by_target_1[class_id],
            by_target_2[class_id],
            junction=(
                junction_by_class[
                    class_id
                ]
            ),
            insert_min=insert_min,
            insert_max=insert_max,
        )

        if geometry is not None:
            class_best[
                class_id
            ] = geometry

    decision = strict_best_class({
        class_id: geometry.score
        for class_id, geometry
        in class_best.items()
    })

    if (
        decision is None
        or decision.supported_class_id
        is None
    ):
        return (
            decision,
            None,
        )

    return (
        decision,
        class_best[
            decision.supported_class_id
        ],
    )


def support_counts(
    decisions: Iterable[
        SupportDecision | None
    ],
) -> Counter[str]:
    """
    Count only strict unique support decisions.
    """
    result: Counter[str] = Counter()

    for decision in decisions:
        if (
            decision is not None
            and decision.supported_class_id
            is not None
        ):
            result[
                decision.supported_class_id
            ] += 1

    return result


def grouped_sam_records(
    lines: Iterable[str],
):
    """
    Yield consecutive SAM records grouped by query name.

    minimap2 emits all alignments for a query together for
    the BGCPhaser short-read invocation. Headers are skipped.
    """
    current_name: str | None = None
    records: list[Alignment] = []

    for raw in lines:
        if raw.startswith("@"):
            continue

        record = parse_sam_record(
            raw
        )

        if current_name is None:
            current_name = (
                record.qname
            )

        if record.qname != current_name:
            yield (
                current_name,
                records,
            )

            current_name = (
                record.qname
            )

            records = []

        records.append(
            record
        )

    if current_name is not None:
        yield (
            current_name,
            records,
        )
