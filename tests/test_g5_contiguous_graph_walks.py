from __future__ import annotations

import importlib.util
from pathlib import Path

import networkx as nx


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarking/scripts/10_localize_g5_contiguous_graph_walks.py"
spec = importlib.util.spec_from_file_location("g5_a1_1", SCRIPT)
assert spec is not None and spec.loader is not None
g5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g5)


def add_state(graph: nx.DiGraph, segment: str, orientation: str, sequence: str) -> tuple[str, str]:
    state = (segment, orientation)
    graph.add_node(state, sequence=sequence)
    return state


def add_edge(graph: nx.DiGraph, left, right, overlap: str = "2M") -> None:
    graph.add_edge(left, right, overlap=overlap, overlaps=[overlap])


def test_split_record_at_missing_edges_preserves_maximal_contiguous_chunks() -> None:
    graph = nx.DiGraph()
    a = add_state(graph, "A", "+", "AAAA")
    b = add_state(graph, "B", "+", "AACC")
    c = add_state(graph, "C", "+", "CCGG")
    d = add_state(graph, "D", "+", "GGTT")
    add_edge(graph, a, b)
    add_edge(graph, c, d)

    chunks, missing = g5.split_record_at_missing_edges(graph, [a, b, c, d])

    assert chunks == [[a, b], [c, d]]
    assert missing == [(b, c)]


def test_reverse_chunk_alignment_collapses_to_same_physical_localization_key() -> None:
    graph = nx.DiGraph()
    a_plus = add_state(graph, "A", "+", "AACCGG")
    b_plus = add_state(graph, "B", "+", "GGTTAA")
    b_minus = add_state(graph, "B", "-", "TTAACC")
    a_minus = add_state(graph, "A", "-", "CCGGTT")
    add_edge(graph, a_plus, b_plus)
    add_edge(graph, b_minus, a_minus)

    plus_path = [a_plus, b_plus]
    reverse_equivalent_path = [b_minus, a_minus]

    plus_row = {
        "tlen": 10,
        "strand": "+",
        "tstart": 1,
        "tend": 9,
    }
    reverse_row = {
        "tlen": 10,
        "strand": "-",
        "tstart": 1,
        "tend": 9,
    }

    plus_loc = g5.localize_alignment(graph, plus_path, plus_row)
    reverse_loc = g5.localize_alignment(graph, reverse_equivalent_path, reverse_row)

    assert plus_loc["key"] == reverse_loc["key"]
    assert plus_loc["subpath"] == [a_plus, b_plus]
    assert reverse_loc["subpath"] == [a_plus, b_plus]
    assert plus_loc["start_offset"] == reverse_loc["start_offset"]
    assert plus_loc["end_offset"] == reverse_loc["end_offset"]


def test_same_inner_state_coordinate_order_requires_left_before_right() -> None:
    state = ("S", "+")
    left = {
        "outer_end_state": state,
        "end_offset": 250,
    }
    right = {
        "outer_start_state": state,
        "start_offset": 500,
    }
    compatible, status = g5.same_inner_state_coordinate_order(left, right)
    assert compatible is True
    assert status == "PASS"

    right_bad = {
        "outer_start_state": state,
        "start_offset": 200,
    }
    compatible, status = g5.same_inner_state_coordinate_order(left, right_bad)
    assert compatible is False
    assert status == "FAIL"

    other_state = ("T", "+")
    right_other = {
        "outer_start_state": other_state,
        "start_offset": 0,
    }
    compatible, status = g5.same_inner_state_coordinate_order(left, right_other)
    assert compatible is True
    assert status == "NOT_APPLICABLE"


def test_primary_contribution_depends_on_role_and_g5_status() -> None:
    assert (
        g5.primary_contribution("PRIMARY_VALIDATION", "PASS")
        == "PENDING_LATER_GATES"
    )
    assert g5.primary_contribution("PRIMARY_VALIDATION", "FAIL") == "NO_G5_FAILURE"
    assert (
        g5.primary_contribution("DEVELOPMENT_SENTINEL_PROTOCOL_AMENDMENT_A1_1", "PASS")
        == "NO_PROTOCOL_AMENDMENT_SENTINEL"
    )
