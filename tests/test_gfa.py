from __future__ import annotations

from pathlib import Path

import pytest

from bgcphaser.gfa import (
    enumerate_oriented_paths,
    enumerate_paths,
    flip_orientation,
    format_oriented_path,
    oriented_node,
    read_gfa,
    reverse_complement,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_reverse_complement_and_orientation_flip() -> None:
    assert reverse_complement("ACGTRYSWKMBDHVN") == "NBDHVKMWSRYACGT"
    assert flip_orientation("+") == "-"
    assert flip_orientation("-") == "+"
    with pytest.raises(ValueError):
        flip_orientation("?")


def test_read_gfa_expands_bidirected_orientation_states(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "orientation.gfa",
        "H\tVN:Z:1.0\nS\tA\tAGTC\nS\tB\tCCAA\nL\tA\t+\tB\t-\t2M\n",
    )
    graph = read_gfa(path)
    assert graph.graph["orientation_model"] == "BIDIRECTED_ORIENTED_STATE"
    assert graph.graph["gfa_segment_count"] == 2
    assert graph.graph["oriented_node_count"] == 4
    assert graph.graph["gfa_link_record_count"] == 1
    assert graph.graph["gfa_non_plus_plus_link_record_count"] == 1
    assert graph.nodes[("A", "+")]["sequence"] == "AGTC"
    assert graph.nodes[("A", "-")]["sequence"] == "GACT"
    assert graph.nodes[("B", "+")]["sequence"] == "CCAA"
    assert graph.nodes[("B", "-")]["sequence"] == "TTGG"
    assert graph.has_edge(("A", "+"), ("B", "-"))
    assert graph.has_edge(("B", "+"), ("A", "-"))


def test_enumerate_two_mixed_orientation_paths(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "two_paths.gfa",
        (
            "S\tS\tAAAA\nS\tA\tCCCC\nS\tB\tGGGG\nS\tT\tTTTT\n"
            "L\tS\t+\tA\t+\t0M\nL\tA\t+\tT\t+\t0M\n"
            "L\tS\t+\tB\t-\t0M\nL\tB\t-\tT\t+\t0M\n"
        ),
    )
    graph = read_gfa(path)
    observed = {
        tuple(candidate)
        for candidate in enumerate_oriented_paths(
            graph,
            ("S", "+"),
            ("T", "+"),
        )
    }
    expected = {
        (("S", "+"), ("A", "+"), ("T", "+")),
        (("S", "+"), ("B", "-"), ("T", "+")),
    }
    assert observed == expected


def test_legacy_plus_orientation_api_is_preserved(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "legacy.gfa",
        "S\tS\tAAAA\nS\tA\tCCCC\nS\tT\tTTTT\nL\tS\t+\tA\t+\t0M\nL\tA\t+\tT\t+\t0M\n",
    )
    assert enumerate_paths(read_gfa(path), "S", "T") == [["S", "A", "T"]]


def test_same_oriented_anchor_returns_trivial_path(tmp_path: Path) -> None:
    path = _write(tmp_path / "same_anchor.gfa", "S\tedge_1\tAACCGGTT\n")
    paths = enumerate_oriented_paths(
        read_gfa(path),
        ("edge_1", "-"),
        ("edge_1", "-"),
    )
    assert paths == [[("edge_1", "-")]]
    assert format_oriented_path(paths[0]) == "edge_1-"


def test_link_to_missing_segment_is_rejected(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "missing.gfa",
        "S\tA\tAAAA\nL\tA\t+\tMISSING\t+\t0M\n",
    )
    with pytest.raises(ValueError, match="missing target segment"):
        read_gfa(path)


def test_oriented_node_rejects_invalid_orientation() -> None:
    with pytest.raises(ValueError):
        oriented_node("A", "?")
