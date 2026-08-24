from bgcphaser.chemistry import chemistry_from_counts
from bgcphaser.scoring import format_q12, rank_candidates, score_candidate
from bgcphaser.sequence import sequence_composite


def test_frozen_bgrr0055_feature_arithmetic():
    sequence = sequence_composite(
        R="0.804347826087",
        L="0.217391304348",
        V="0.860371272914",
        H="0.043478260870",
    )
    assert format_q12(sequence.score) == "0.709658035620"

    chemistry = chemistry_from_counts(
        c_supported=30,
        c_generic_unresolved=1,
        c_prediction_missing=0,
        c_assessable=31,
        m_complete=36,
        m_incomplete=1,
        m_total=37,
    )
    assert format_q12(chemistry.C) == "0.967741935484"
    assert format_q12(chemistry.M) == "0.972972972973"

    combined = score_candidate(
        "BGRR0055_walk_00001",
        sequence.score,
        chemistry.C,
        chemistry.M,
    )
    assert format_q12(combined.chemistry_score) == "0.970357454228"
    assert format_q12(combined.combined_score) == "0.840007744924"


def test_frozen_bgrr0074_control_score():
    control = score_candidate(
        "BGRR0074_walk_00001",
        "0.735152428026",
        "1.000000000000",
        "0.812500000000",
    )
    assert format_q12(control.chemistry_score) == "0.906250000000"
    assert format_q12(control.combined_score) == "0.820701214013"

    ranked = rank_candidates([control])
    assert ranked[0].rank == 1
    assert ranked[0].status == "RANKED"


def test_missing_required_component_remains_unranked():
    missing = rank_candidates(
        [
            score_candidate(
                "missing_C",
                "0.800000000000",
                "NA",
                "0.900000000000",
            )
        ]
    )[0]
    assert missing.rank is None
    assert missing.combined_score is None
    assert missing.status == "UNRANKED_C_NA"
