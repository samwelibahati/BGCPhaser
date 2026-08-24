from __future__ import annotations

from dataclasses import dataclass, replace
from decimal import (
    Decimal,
    InvalidOperation,
    ROUND_HALF_EVEN,
)
from typing import Iterable


Q12 = Decimal("0.000000000001")
HALF = Decimal("0.5")

NA_TOKENS = {
    "",
    "NA",
    "N/A",
    "NONE",
    "NULL",
    ".",
}


def q12(value: Decimal) -> Decimal:
    return value.quantize(
        Q12,
        rounding=ROUND_HALF_EVEN,
    )


def parse_score(
    value: str | Decimal | None,
    *,
    name: str,
    allow_na: bool = False,
) -> Decimal | None:
    if value is None:
        if allow_na:
            return None
        raise ValueError(f"{name} is required")

    text = str(value).strip()

    if text.upper() in NA_TOKENS:
        if allow_na:
            return None
        raise ValueError(f"{name} is required")

    try:
        result = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError(
            f"{name} is not a valid decimal: {text}"
        ) from exc

    if not result.is_finite():
        raise ValueError(
            f"{name} must be finite"
        )

    if not (
        Decimal("0")
        <= result
        <= Decimal("1")
    ):
        raise ValueError(
            f"{name} must be between 0 and 1"
        )

    return result


def format_q12(
    value: Decimal | None,
) -> str:
    if value is None:
        return "NA"

    return f"{value:.12f}"


@dataclass(frozen=True)
class CandidateScore:
    candidate_id: str
    sequence_score: Decimal | None
    chemistry_c: Decimal | None
    chemistry_m: Decimal | None
    chemistry_score: Decimal | None
    combined_score: Decimal | None
    rank: int | None = None
    status: str = "UNRANKED"


def score_candidate(
    candidate_id: str,
    sequence_score: str | Decimal | None,
    chemistry_c: str | Decimal | None,
    chemistry_m: str | Decimal | None,
) -> CandidateScore:
    """
    Frozen complete-case combined scoring.

    S_chem =
        q12(0.5*C + 0.5*M)

    S_combined =
        q12(0.5*S_seq + 0.5*S_chem)

    Any required NA yields S_combined=NA.
    No missing weight is redistributed.
    """
    candidate_id = str(
        candidate_id
    ).strip()

    if not candidate_id:
        raise ValueError(
            "candidate_id is required"
        )

    seq = parse_score(
        sequence_score,
        name="sequence_score",
        allow_na=True,
    )

    c_value = parse_score(
        chemistry_c,
        name="chemistry_c",
        allow_na=True,
    )

    m_value = parse_score(
        chemistry_m,
        name="chemistry_m",
        allow_na=True,
    )

    if (
        c_value is None
        or m_value is None
    ):
        chemistry_score = None
    else:
        chemistry_score = q12(
            HALF * c_value
            + HALF * m_value
        )

    if (
        seq is None
        or chemistry_score is None
    ):
        combined_score = None
    else:
        combined_score = q12(
            HALF * seq
            + HALF * chemistry_score
        )

    return CandidateScore(
        candidate_id=candidate_id,
        sequence_score=seq,
        chemistry_c=c_value,
        chemistry_m=m_value,
        chemistry_score=chemistry_score,
        combined_score=combined_score,
    )


def _unranked_status(
    item: CandidateScore,
) -> str:
    if item.sequence_score is None:
        return (
            "UNRANKED_SEQUENCE_SCORE_NA"
        )

    if (
        item.chemistry_c is None
        and item.chemistry_m is None
    ):
        return "UNRANKED_C_AND_M_NA"

    if item.chemistry_c is None:
        return "UNRANKED_C_NA"

    if item.chemistry_m is None:
        return "UNRANKED_M_NA"

    return (
        "UNRANKED_REQUIRED_COMPONENT_NA"
    )


def rank_candidates(
    candidates: Iterable[CandidateScore],
) -> list[CandidateScore]:
    """
    Standard competition ranking.

    Exact q12 ties share a rank.
    Candidate ID never breaks a scientific tie.
    """
    values = list(candidates)

    ids = [
        item.candidate_id
        for item in values
    ]

    if len(ids) != len(set(ids)):
        raise ValueError(
            "duplicate candidate_id"
        )

    counts: dict[Decimal, int] = {}

    for item in values:
        if item.combined_score is None:
            continue

        counts[item.combined_score] = (
            counts.get(
                item.combined_score,
                0,
            )
            + 1
        )

    rank_by_score: dict[
        Decimal,
        int,
    ] = {}

    next_rank = 1

    for value in sorted(
        counts,
        reverse=True,
    ):
        rank_by_score[value] = next_rank
        next_rank += counts[value]

    ranked = []

    for item in values:
        if item.combined_score is None:
            ranked.append(
                replace(
                    item,
                    rank=None,
                    status=_unranked_status(
                        item
                    ),
                )
            )
        else:
            ranked.append(
                replace(
                    item,
                    rank=rank_by_score[
                        item.combined_score
                    ],
                    status="RANKED",
                )
            )

    ranked.sort(
        key=lambda item: (
            item.rank is None,
            (
                item.rank
                if item.rank is not None
                else 10**18
            ),
            item.candidate_id,
        )
    )

    return ranked
