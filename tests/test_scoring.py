from decimal import Decimal

from bgcphaser.scoring import (
    format_q12,
    q12,
    rank_candidates,
    score_candidate,
)


def test_q12_uses_round_half_even():
    assert q12(Decimal("0.1234567890125")) == Decimal("0.123456789012")
    assert q12(Decimal("0.1234567890135")) == Decimal("0.123456789014")


def test_frozen_staged_scoring():
    result = score_candidate(
        "candidate_A",
        "0.600000000000",
        "0.800000000000",
        "1.000000000000",
    )
    assert format_q12(result.chemistry_score) == "0.900000000000"
    assert format_q12(result.combined_score) == "0.750000000000"


def test_competition_ranking_and_exact_ties():
    candidates = [
        score_candidate(
            "candidate_A",
            "0.600000000000",
            "0.800000000000",
            "1.000000000000",
        ),
        score_candidate(
            "candidate_B",
            "0.700000000000",
            "0.800000000000",
            "0.800000000000",
        ),
        score_candidate(
            "candidate_C",
            "0.500000000000",
            "0.500000000000",
            "0.500000000000",
        ),
        score_candidate(
            "candidate_D",
            "NA",
            "0.900000000000",
            "0.900000000000",
        ),
    ]
    ranked = {
        item.candidate_id: item
        for item in rank_candidates(candidates)
    }
    assert ranked["candidate_A"].rank == 1
    assert ranked["candidate_B"].rank == 1
    assert ranked["candidate_C"].rank == 3
    assert ranked["candidate_D"].rank is None
    assert ranked["candidate_D"].status == "UNRANKED_SEQUENCE_SCORE_NA"


def test_missing_chemistry_is_not_redistributed():
    result = score_candidate(
        "candidate_A",
        "0.700000000000",
        "NA",
        "0.900000000000",
    )
    assert result.chemistry_score is None
    assert result.combined_score is None
    ranked = rank_candidates([result])[0]
    assert ranked.rank is None
    assert ranked.status == "UNRANKED_C_NA"
