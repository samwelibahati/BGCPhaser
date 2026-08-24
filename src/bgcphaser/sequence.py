from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
import statistics
from typing import Mapping, Sequence

from .scoring import parse_score, q12


OrientedNode = tuple[str, str]
OrientedEdge = tuple[OrientedNode, OrientedNode]


@dataclass(frozen=True)
class SequenceComposite:
    R: Decimal | None
    L: Decimal | None
    V: Decimal | None
    H: Decimal | None
    one_minus_H: Decimal | None
    defined_components: tuple[str, ...]
    score: Decimal | None


@dataclass(frozen=True)
class CoverageConsistency:
    V: Decimal | None
    mad: Decimal | None
    usable_unique_segment_count: int


@dataclass(frozen=True)
class AssemblerDiscordance:
    H: Decimal | None
    discordant_transition_count: int | None
    transition_denominator: int


def support_fraction(
    supported: int,
    assessable: int,
) -> Decimal | None:
    """
    Frozen R/L candidate aggregation.

    R or L =
        supported assessable transition occurrences
        / assessable transition occurrences.
    """
    if assessable < 0:
        raise ValueError("assessable must be >= 0")

    if supported < 0 or supported > assessable:
        raise ValueError(
            "supported must be between 0 and assessable"
        )

    if assessable == 0:
        return None

    with localcontext() as context:
        context.prec = 50
        value = (
            Decimal(supported)
            / Decimal(assessable)
        )

    return q12(value)


def coverage_consistency(
    walk: Sequence[OrientedNode],
    depth_by_segment: Mapping[str, float],
) -> CoverageConsistency:
    """
    Frozen V definition.

    For each unique segment s:
        x_s = log2(depth_s / multiplicity_s)

    V = 1 / (1 + MAD(x_s))

    Orientation is collapsed for depth.
    """
    multiplicity = Counter(
        segment_id
        for segment_id, _ in walk
    )

    copy_adjusted_logs: list[float] = []

    for segment_id in sorted(multiplicity):
        if segment_id not in depth_by_segment:
            continue

        depth = float(
            depth_by_segment[segment_id]
        )

        per_copy_depth = (
            depth
            / multiplicity[segment_id]
        )

        if per_copy_depth <= 0:
            continue

        copy_adjusted_logs.append(
            math.log2(per_copy_depth)
        )

    usable = len(copy_adjusted_logs)

    if usable < 2:
        return CoverageConsistency(
            V=None,
            mad=None,
            usable_unique_segment_count=usable,
        )

    center = statistics.median(
        copy_adjusted_logs
    )

    mad_float = statistics.median(
        [
            abs(value - center)
            for value in copy_adjusted_logs
        ]
    )

    v_float = 1.0 / (1.0 + mad_float)

    return CoverageConsistency(
        V=Decimal(f"{v_float:.12f}"),
        mad=Decimal(f"{mad_float:.12f}"),
        usable_unique_segment_count=usable,
    )


def assembler_discordance(
    walk: Sequence[OrientedNode],
    assembler_adjacency: set[OrientedEdge] | None,
) -> AssemblerDiscordance:
    """
    Frozen H definition.

    H =
        discordant oriented transition occurrences
        / all candidate transition occurrences.
    """
    transitions = list(
        zip(
            walk,
            walk[1:],
        )
    )

    denominator = len(transitions)

    if (
        denominator == 0
        or not assembler_adjacency
    ):
        return AssemblerDiscordance(
            H=None,
            discordant_transition_count=None,
            transition_denominator=denominator,
        )

    discordant = sum(
        (left, right)
        not in assembler_adjacency
        for left, right in transitions
    )

    h_float = (
        discordant
        / denominator
    )

    return AssemblerDiscordance(
        H=Decimal(f"{h_float:.12f}"),
        discordant_transition_count=discordant,
        transition_denominator=denominator,
    )


def sequence_composite(
    R: str | Decimal | None,
    L: str | Decimal | None,
    V: str | Decimal | None,
    H: str | Decimal | None,
) -> SequenceComposite:
    """
    Frozen sequence-score definition.

    S_seq =
        arithmetic mean of all defined components
        among R, L, V, and 1-H.

    Missing components are excluded from the mean.
    They are never imputed as zero.
    """
    r = parse_score(
        R,
        name="R",
        allow_na=True,
    )

    l = parse_score(
        L,
        name="L",
        allow_na=True,
    )

    v = parse_score(
        V,
        name="V",
        allow_na=True,
    )

    h = parse_score(
        H,
        name="H",
        allow_na=True,
    )

    if (
        v is not None
        and v <= Decimal("0")
    ):
        raise ValueError(
            "defined V must be > 0"
        )

    one_minus_h = (
        q12(
            Decimal("1") - h
        )
        if h is not None
        else None
    )

    labelled = (
        ("R", r),
        ("L", l),
        ("V", v),
        ("1-H", one_minus_h),
    )

    components = [
        value
        for _, value in labelled
        if value is not None
    ]

    labels = tuple(
        label
        for label, value in labelled
        if value is not None
    )

    if not components:
        score = None
    else:
        with localcontext() as context:
            context.prec = 50
            value = (
                sum(
                    components,
                    Decimal("0"),
                )
                / Decimal(len(components))
            )

        score = q12(value)

    return SequenceComposite(
        R=r,
        L=l,
        V=v,
        H=h,
        one_minus_H=one_minus_h,
        defined_components=labels,
        score=score,
    )
