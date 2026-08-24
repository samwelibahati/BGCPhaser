from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal


@dataclass(frozen=True)
class ChemistryFeatures:
    C: Decimal | None
    C_status: str
    C_supported_module_count: int
    C_generic_unresolved_module_count: int
    C_prediction_missing_module_count: int
    C_assessable_explicit_loader_module_count: int

    M: Decimal | None
    M_status: str
    M_complete_module_count: int
    M_incomplete_module_count: int
    M_total_module_count: int


def _fraction12(
    numerator: int,
    denominator: int,
) -> Decimal:
    """
    Preserve the locked chemistry extraction
    arithmetic: Python ratio formatted to 12 places.
    """
    return Decimal(
        f"{numerator / denominator:.12f}"
    )


def chemistry_from_counts(
    *,
    c_supported: int,
    c_generic_unresolved: int,
    c_prediction_missing: int,
    c_assessable: int,
    m_complete: int,
    m_incomplete: int,
    m_total: int,
) -> ChemistryFeatures:
    """
    Frozen candidate-level C and M aggregation.

    C =
        supported assessable explicit-loader modules
        / all assessable explicit-loader modules.

    M =
        complete NRPS/PKS modules
        / all NRPS/PKS module objects.
    """
    counts = (
        c_supported,
        c_generic_unresolved,
        c_prediction_missing,
        c_assessable,
        m_complete,
        m_incomplete,
        m_total,
    )

    if any(value < 0 for value in counts):
        raise ValueError(
            "chemistry counts must be >= 0"
        )

    if (
        c_supported
        + c_generic_unresolved
        + c_prediction_missing
        != c_assessable
    ):
        raise ValueError(
            "C evidence counts do not sum "
            "to assessable module count"
        )

    if (
        m_complete
        + m_incomplete
        != m_total
    ):
        raise ValueError(
            "M complete/incomplete counts "
            "do not sum to total module count"
        )

    if c_assessable == 0:
        c_value = None
        c_status = (
            "NA_ZERO_ASSESSABLE_"
            "EXPLICIT_LOADER_MODULES"
        )

    elif c_prediction_missing > 0:
        c_value = None
        c_status = (
            "NA_ASSESSABLE_MODULE_"
            "PREDICTION_MISSING"
        )

    else:
        c_value = _fraction12(
            c_supported,
            c_assessable,
        )
        c_status = "DEFINED"

    if m_total == 0:
        m_value = None
        m_status = (
            "NA_ZERO_NRPS_PKS_MODULES"
        )

    else:
        m_value = _fraction12(
            m_complete,
            m_total,
        )
        m_status = "DEFINED"

    return ChemistryFeatures(
        C=c_value,
        C_status=c_status,
        C_supported_module_count=c_supported,
        C_generic_unresolved_module_count=(
            c_generic_unresolved
        ),
        C_prediction_missing_module_count=(
            c_prediction_missing
        ),
        C_assessable_explicit_loader_module_count=(
            c_assessable
        ),
        M=m_value,
        M_status=m_status,
        M_complete_module_count=m_complete,
        M_incomplete_module_count=m_incomplete,
        M_total_module_count=m_total,
    )


LOADER_TYPES = {
    "AMP-binding",
    "PKS_AT",
}

GENERIC_UNRESOLVED = {
    "",
    "X",
    "(UNKNOWN)",
    "UNKNOWN",
    "NRP",
    "PK",
}


def _string_values(
    value,
    *,
    context: str,
) -> list[str]:
    if value is None:
        return []

    if isinstance(
        value,
        str,
    ):
        return [value]

    if isinstance(
        value,
        list,
    ):
        if not all(
            isinstance(item, str)
            for item in value
        ):
            raise ValueError(
                f"{context}: expected list[str]"
            )

        return list(value)

    raise ValueError(
        f"{context}: unexpected "
        f"value type "
        f"{type(value).__name__}"
    )


def _is_generic_consensus(
    value: str,
) -> bool:
    return (
        value.strip().upper()
        in GENERIC_UNRESOLVED
    )


def chemistry_from_antismash_data(
    data: dict,
) -> ChemistryFeatures:
    """
    Extract frozen BGCPhaser C/M evidence from an
    antiSMASH JSON object.

    This consumes only intrinsic NRPS/PKS annotation
    evidence. ClusterBlast/ClusterCompare outputs are
    explicitly rejected as scoring inputs.
    """
    if not isinstance(
        data,
        dict,
    ):
        raise ValueError(
            "antiSMASH JSON top level "
            "must be an object"
        )

    records = data.get(
        "records"
    )

    if (
        not isinstance(records, list)
        or len(records) != 1
    ):
        raise ValueError(
            "expected exactly one "
            "antiSMASH record"
        )

    record = records[0]

    if not isinstance(
        record,
        dict,
    ):
        raise ValueError(
            "antiSMASH record "
            "must be an object"
        )

    record_modules = record.get(
        "modules"
    )

    if not isinstance(
        record_modules,
        dict,
    ):
        raise ValueError(
            "serialized record.modules "
            "is absent"
        )

    for prohibited in (
        "antismash.modules.clusterblast",
        "antismash.modules.cluster_compare",
    ):
        if prohibited in record_modules:
            raise ValueError(
                "reference-comparison result "
                "cannot be a scoring input: "
                f"{prohibited}"
            )

    nrps_results = record_modules.get(
        "antismash.modules.nrps_pks"
    )

    if not isinstance(
        nrps_results,
        dict,
    ):
        raise ValueError(
            "NRPS/PKS result object absent"
        )

    consensus = nrps_results.get(
        "consensus"
    )

    if not isinstance(
        consensus,
        dict,
    ):
        raise ValueError(
            "NRPS/PKS consensus "
            "must be a dict"
        )

    if not all(
        isinstance(key, str)
        and isinstance(value, str)
        for key, value
        in consensus.items()
    ):
        raise ValueError(
            "NRPS/PKS consensus schema drift"
        )

    features = record.get(
        "features"
    )

    if not isinstance(
        features,
        list,
    ):
        raise ValueError(
            "serialized antiSMASH "
            "features absent"
        )

    domain_index: dict[
        str,
        list[str],
    ] = {}

    for feature in features:
        if (
            not isinstance(feature, dict)
            or feature.get("type")
            != "aSDomain"
        ):
            continue

        qualifiers = feature.get(
            "qualifiers"
        )

        if not isinstance(
            qualifiers,
            dict,
        ):
            raise ValueError(
                "aSDomain qualifiers absent"
            )

        domain_ids = _string_values(
            qualifiers.get(
                "domain_id"
            ),
            context="aSDomain.domain_id",
        )

        if len(domain_ids) != 1:
            raise ValueError(
                "aSDomain must have exactly "
                "one domain_id"
            )

        domain_types = _string_values(
            qualifiers.get(
                "aSDomain"
            ),
            context="aSDomain.aSDomain",
        )

        if not domain_types:
            raise ValueError(
                "aSDomain lacks domain type"
            )

        domain_id = domain_ids[0]

        if domain_id in domain_index:
            raise ValueError(
                "duplicate aSDomain ID: "
                f"{domain_id}"
            )

        domain_index[
            domain_id
        ] = domain_types

    modules = [
        feature
        for feature in features
        if (
            isinstance(feature, dict)
            and feature.get("type")
            == "aSModule"
        )
    ]

    c_assessable = 0
    c_supported = 0
    c_generic = 0
    c_missing = 0

    m_complete = 0
    m_incomplete = 0

    for module_index, feature in enumerate(
        modules,
        start=1,
    ):
        qualifiers = feature.get(
            "qualifiers"
        )

        if not isinstance(
            qualifiers,
            dict,
        ):
            raise ValueError(
                "aSModule qualifiers absent"
            )

        has_complete = (
            "complete" in qualifiers
        )

        has_incomplete = (
            "incomplete" in qualifiers
        )

        if (
            has_complete
            == has_incomplete
        ):
            raise ValueError(
                f"module {module_index} "
                "must carry exactly one of "
                "complete/incomplete"
            )

        complete = has_complete

        if complete:
            m_complete += 1
        else:
            m_incomplete += 1

        domain_ids = _string_values(
            qualifiers.get(
                "domains"
            ),
            context=(
                f"module "
                f"{module_index}.domains"
            ),
        )

        loader_entries: list[
            tuple[str, str]
        ] = []

        for domain_id in domain_ids:
            if domain_id not in domain_index:
                raise ValueError(
                    f"module {module_index} "
                    "references unknown domain "
                    f"{domain_id}"
                )

            for domain_type in (
                domain_index[
                    domain_id
                ]
            ):
                if (
                    domain_type
                    in LOADER_TYPES
                ):
                    loader_entries.append(
                        (
                            domain_id,
                            domain_type,
                        )
                    )

        assessable = (
            complete
            and bool(loader_entries)
        )

        if not assessable:
            continue

        c_assessable += 1

        loader_statuses: list[
            str
        ] = []

        for (
            domain_id,
            _domain_type,
        ) in loader_entries:

            if domain_id not in consensus:
                loader_statuses.append(
                    "MISSING_CONSENSUS"
                )

            elif _is_generic_consensus(
                consensus[domain_id]
            ):
                loader_statuses.append(
                    "GENERIC_UNRESOLVED"
                )

            else:
                loader_statuses.append(
                    "NONGENERIC_SUPPORTED"
                )

        if (
            "NONGENERIC_SUPPORTED"
            in loader_statuses
        ):
            c_supported += 1

        elif (
            "MISSING_CONSENSUS"
            in loader_statuses
        ):
            c_missing += 1

        else:
            c_generic += 1

    return chemistry_from_counts(
        c_supported=c_supported,
        c_generic_unresolved=c_generic,
        c_prediction_missing=c_missing,
        c_assessable=c_assessable,
        m_complete=m_complete,
        m_incomplete=m_incomplete,
        m_total=len(modules),
    )


def chemistry_from_antismash_json(
    path,
) -> ChemistryFeatures:
    """
    Load one antiSMASH JSON file and extract C/M.
    """
    import json
    from pathlib import Path

    json_path = Path(path)

    with json_path.open(
        encoding="utf-8"
    ) as handle:
        data = json.load(handle)

    return chemistry_from_antismash_data(
        data
    )
