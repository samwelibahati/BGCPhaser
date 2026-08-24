from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

import networkx as nx


OrientedNode: TypeAlias = tuple[str, str]
OrientedPath: TypeAlias = list[OrientedNode]


_RC_TABLE = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return sequence.translate(_RC_TABLE)[::-1]


def flip_orientation(orientation: str) -> str:
    """Return the opposite GFA segment orientation."""
    if orientation == "+":
        return "-"
    if orientation == "-":
        return "+"

    raise ValueError(
        f"Invalid GFA orientation: {orientation!r}"
    )


def oriented_node(
    segment_id: str,
    orientation: str,
) -> OrientedNode:
    """Construct and validate an oriented graph node."""
    if orientation not in {"+", "-"}:
        raise ValueError(
            f"Invalid GFA orientation: {orientation!r}"
        )

    return str(segment_id), orientation


def format_oriented_node(
    node: OrientedNode,
) -> str:
    """Format an oriented node for human-readable output."""
    segment_id, orientation = node
    return f"{segment_id}{orientation}"


def format_oriented_path(
    path: OrientedPath,
) -> str:
    """Format an oriented path for audit output."""
    return "->".join(
        format_oriented_node(node)
        for node in path
    )


def path_segment_ids(
    path: OrientedPath,
) -> list[str]:
    """Project an oriented path onto segment identifiers."""
    return [
        segment_id
        for segment_id, _orientation in path
    ]


def _add_transition(
    graph: nx.DiGraph,
    source: OrientedNode,
    target: OrientedNode,
    overlap: str,
    provenance: str,
) -> None:
    """Add one directed oriented transition without duplicating topology."""
    if graph.has_edge(
        source,
        target,
    ):
        attrs = graph.edges[
            source,
            target,
        ]

        overlaps = list(
            attrs.get(
                "overlaps",
                (),
            )
        )

        if overlap not in overlaps:
            overlaps.append(
                overlap
            )

        provenances = list(
            attrs.get(
                "provenances",
                (),
            )
        )

        if provenance not in provenances:
            provenances.append(
                provenance
            )

        attrs[
            "overlaps"
        ] = tuple(
            overlaps
        )

        attrs[
            "provenances"
        ] = tuple(
            provenances
        )

        attrs[
            "link_observations"
        ] = (
            int(
                attrs.get(
                    "link_observations",
                    1,
                )
            )
            + 1
        )

        return

    graph.add_edge(
        source,
        target,
        overlap=overlap,
        overlaps=(overlap,),
        provenances=(provenance,),
        link_observations=1,
    )


def read_gfa(
    path: str | Path,
) -> nx.DiGraph:
    """Read GFA1 S/L topology as a bidirected oriented-state graph.

    Every GFA segment is represented by two graph nodes and each GFA link
    contributes the stated traversal plus its reverse-complement reciprocal.
    """
    gfa_path = Path(path)

    segment_records: dict[
        str,
        tuple[str | None, tuple[str, ...]],
    ] = {}
    link_records: list[
        tuple[str, str, str, str, str, int]
    ] = []
    gfa_path_records = 0
    non_plus_plus_links = 0

    with gfa_path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            record_type = fields[0]

            if record_type == "S":
                if len(fields) < 3:
                    raise ValueError(
                        f"Malformed GFA S record at {gfa_path}:{line_number}"
                    )
                segment_id = fields[1]
                raw_sequence = fields[2]
                if segment_id in segment_records:
                    raise ValueError(
                        f"Duplicate GFA segment ID {segment_id!r} at {gfa_path}:{line_number}"
                    )
                sequence = None if raw_sequence == "*" else raw_sequence
                segment_records[segment_id] = (sequence, tuple(fields[3:]))

            elif record_type == "L":
                if len(fields) < 6:
                    raise ValueError(
                        f"Malformed GFA L record at {gfa_path}:{line_number}"
                    )
                source, source_orientation = fields[1], fields[2]
                target, target_orientation = fields[3], fields[4]
                overlap = fields[5]
                oriented_node(source, source_orientation)
                oriented_node(target, target_orientation)
                if source_orientation != "+" or target_orientation != "+":
                    non_plus_plus_links += 1
                link_records.append(
                    (source, source_orientation, target, target_orientation, overlap, line_number)
                )
            elif record_type in {"P", "W"}:
                gfa_path_records += 1

    if not segment_records:
        raise ValueError(f"No segments found in GFA: {gfa_path}")

    graph = nx.DiGraph()
    for segment_id, (sequence, raw_tags) in segment_records.items():
        plus = oriented_node(segment_id, "+")
        minus = oriented_node(segment_id, "-")
        graph.add_node(
            plus,
            segment_id=segment_id,
            orientation="+",
            sequence=sequence,
            raw_tags=raw_tags,
        )
        graph.add_node(
            minus,
            segment_id=segment_id,
            orientation="-",
            sequence=None if sequence is None else reverse_complement(sequence),
            raw_tags=raw_tags,
        )

    for source, source_orientation, target, target_orientation, overlap, line_number in link_records:
        if source not in segment_records:
            raise ValueError(
                f"GFA link at {gfa_path}:{line_number} references missing source segment {source!r}"
            )
        if target not in segment_records:
            raise ValueError(
                f"GFA link at {gfa_path}:{line_number} references missing target segment {target!r}"
            )
        stated_source = oriented_node(source, source_orientation)
        stated_target = oriented_node(target, target_orientation)
        reciprocal_source = oriented_node(target, flip_orientation(target_orientation))
        reciprocal_target = oriented_node(source, flip_orientation(source_orientation))
        _add_transition(graph, stated_source, stated_target, overlap, "stated")
        _add_transition(
            graph,
            reciprocal_source,
            reciprocal_target,
            overlap,
            "reverse_complement_reciprocal",
        )

    graph.graph["source_path"] = str(gfa_path)
    graph.graph["orientation_model"] = "BIDIRECTED_ORIENTED_STATE"
    graph.graph["gfa_segment_count"] = len(segment_records)
    graph.graph["oriented_node_count"] = graph.number_of_nodes()
    graph.graph["gfa_link_record_count"] = len(link_records)
    graph.graph["directed_transition_count"] = graph.number_of_edges()
    graph.graph["gfa_path_record_count"] = gfa_path_records
    graph.graph["gfa_non_plus_plus_link_record_count"] = non_plus_plus_links
    return graph


def enumerate_oriented_paths(
    graph: nx.DiGraph,
    start: OrientedNode,
    end: OrientedNode,
    max_paths: int = 10000,
    max_nodes: int = 100,
) -> list[OrientedPath]:
    """Enumerate bounded simple paths between oriented graph states."""
    if max_paths < 1:
        raise ValueError("max_paths must be >= 1")
    if max_nodes < 1:
        raise ValueError("max_nodes must be >= 1")

    start = oriented_node(start[0], start[1])
    end = oriented_node(end[0], end[1])
    if start not in graph or end not in graph:
        raise KeyError(f"Start/end oriented node missing: {start}, {end}")
    if start == end:
        return [[start]]

    paths: list[OrientedPath] = []
    stack: list[tuple[OrientedNode, OrientedPath, set[OrientedNode]]] = [
        (start, [start], {start})
    ]
    while stack:
        node, path, seen = stack.pop()
        if len(path) >= max_nodes:
            continue
        successors = sorted(
            graph.successors(node),
            key=lambda item: (item[0], item[1]),
            reverse=True,
        )
        for next_node in successors:
            if next_node in seen:
                continue
            next_path = [*path, next_node]
            if next_node == end:
                paths.append(next_path)
                if len(paths) >= max_paths:
                    return paths
                continue
            next_seen = set(seen)
            next_seen.add(next_node)
            stack.append((next_node, next_path, next_seen))
    return paths


def enumerate_paths(
    graph: nx.DiGraph,
    start: str,
    end: str,
    max_paths: int = 10000,
    max_nodes: int = 100,
) -> list[list[str]]:
    """Legacy plus-orientation path API used by the prototype scorer."""
    oriented_paths = enumerate_oriented_paths(
        graph,
        oriented_node(start, "+"),
        oriented_node(end, "+"),
        max_paths=max_paths,
        max_nodes=max_nodes,
    )
    return [path_segment_ids(path) for path in oriented_paths]


from collections import Counter as _Counter
from collections.abc import Iterable as _Iterable
from collections.abc import Mapping as _Mapping
from pathlib import Path as _Path

import networkx as _nx


def parse_spades_path_records(
    path: str | _Path,
) -> list[OrientedPath]:
    """Parse SPAdes contigs.paths/scaffolds.paths into oriented subpaths."""
    path = _Path(path)
    records: list[OrientedPath] = []
    current_header: str | None = None
    pieces: list[str] = []

    def flush() -> None:
        if current_header is None:
            return
        text = "".join(piece.strip() for piece in pieces)
        for chunk in text.split(";"):
            chunk = chunk.strip().strip(",")
            if not chunk:
                continue
            oriented: OrientedPath = []
            for token in chunk.split(","):
                token = token.strip()
                if not token:
                    continue
                orientation = token[-1]
                if orientation not in {"+", "-"}:
                    raise ValueError(
                        f"Invalid SPAdes path token in {path}: {token!r}"
                    )
                segment_id = token[:-1]
                if not segment_id:
                    raise ValueError(
                        f"Missing segment ID in {path}: {token!r}"
                    )
                oriented.append((segment_id, orientation))
            if oriented:
                records.append(oriented)

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("NODE_"):
                flush()
                current_header = line
                pieces = []
            else:
                if current_header is None:
                    raise ValueError(
                        "SPAdes path content encountered "
                        f"before NODE header in {path}: {line!r}"
                    )
                pieces.append(line)
    flush()
    return records


def derive_assembler_walk_caps(
    path_files: _Iterable[str | _Path],
) -> tuple[
    dict[OrientedNode, int],
    dict[tuple[OrientedNode, OrientedNode], int],
]:
    """Derive state and directed-edge multiplicity caps from assembler paths."""
    state_caps: dict[OrientedNode, int] = {}
    edge_caps: dict[tuple[OrientedNode, OrientedNode], int] = {}
    for path_file in path_files:
        for assembler_path in parse_spades_path_records(path_file):
            state_counts = _Counter(assembler_path)
            edge_counts = _Counter(zip(assembler_path, assembler_path[1:]))
            for state, count in state_counts.items():
                state_caps[state] = max(state_caps.get(state, 1), count)
            for edge, count in edge_counts.items():
                edge_caps[edge] = max(edge_caps.get(edge, 1), count)
    return state_caps, edge_caps


def derive_anchor_bearing_walk_bound(
    graph,
    path_files: _Iterable[str | _Path],
    start: OrientedNode,
    end: OrientedNode,
) -> int:
    """Derive a reference-independent walk-length bound from assembler paths."""
    if start not in graph:
        raise KeyError(f"Start state not present in graph: {start}")
    if end not in graph:
        raise KeyError(f"End state not present in graph: {end}")
    if start == end:
        return 1

    records: list[OrientedPath] = []
    for path_file in path_files:
        records.extend(parse_spades_path_records(path_file))

    for assembler_path in records:
        for source, target in zip(assembler_path, assembler_path[1:]):
            if not graph.has_edge(source, target):
                raise ValueError(
                    "Assembler path transition not present "
                    f"in selected graph: {source}->{target}"
                )

    left_occurrences: list[tuple[OrientedPath, int]] = []
    right_occurrences: list[tuple[OrientedPath, int]] = []
    best = 0

    for assembler_path in records:
        start_indices = [
            index for index, state in enumerate(assembler_path) if state == start
        ]
        end_indices = [
            index for index, state in enumerate(assembler_path) if state == end
        ]
        for index in start_indices:
            left_occurrences.append((assembler_path, index))
        for index in end_indices:
            right_occurrences.append((assembler_path, index))
        for left_index in start_indices:
            for right_index in end_indices:
                if left_index <= right_index:
                    best = max(best, right_index - left_index + 1)

    if not left_occurrences:
        raise ValueError(
            f"No supplied assembler subpath contains the left anchor state {start}"
        )
    if not right_occurrences:
        raise ValueError(
            f"No supplied assembler subpath contains the right anchor state {end}"
        )

    for left_path, left_index in left_occurrences:
        left_suffix = left_path[left_index:]
        left_exit = left_suffix[-1]
        for right_path, right_index in right_occurrences:
            right_prefix = right_path[:right_index + 1]
            right_entry = right_prefix[0]
            try:
                connector = _nx.shortest_path(graph, left_exit, right_entry)
            except _nx.NetworkXNoPath:
                continue
            if left_exit == right_entry:
                combined_nodes = len(left_suffix) + len(right_prefix) - 1
            else:
                combined_nodes = (
                    len(left_suffix)
                    + max(0, len(connector) - 2)
                    + len(right_prefix)
                )
            best = max(best, combined_nodes)

    if best < 1:
        raise ValueError(
            f"No anchor-bearing assembler template can connect {start} to {end}"
        )
    return best


def enumerate_oriented_walks(
    graph,
    start: OrientedNode,
    end: OrientedNode,
    *,
    state_caps: _Mapping[OrientedNode, int],
    edge_caps: _Mapping[tuple[OrientedNode, OrientedNode], int],
    max_nodes: int,
    max_walks: int = 10_000,
) -> list[OrientedPath]:
    """Enumerate assembler-bounded directed walks between oriented anchors."""
    if max_nodes < 1:
        raise ValueError("max_nodes must be at least 1")
    if max_walks < 1:
        raise ValueError("max_walks must be at least 1")
    if start not in graph:
        raise KeyError(f"Start state not present in graph: {start}")
    if end not in graph:
        raise KeyError(f"End state not present in graph: {end}")
    if start == end:
        return [[start]]

    forward = _nx.descendants(graph, start) | {start}
    backward = _nx.ancestors(graph, end) | {end}
    corridor_nodes = forward & backward
    corridor = graph.subgraph(corridor_nodes)
    if end not in corridor:
        return []

    reverse_corridor = corridor.reverse(copy=False)
    distance_to_end = _nx.single_source_shortest_path_length(
        reverse_corridor,
        end,
    )
    initial_state_visits = _Counter({start: 1})
    stack = [(start, [start], initial_state_visits, _Counter())]
    walks: list[OrientedPath] = []

    while stack:
        node, path, state_visits, edge_visits = stack.pop()
        if len(path) >= max_nodes:
            continue
        successors = sorted(
            corridor.successors(node),
            key=lambda state: (state[0], state[1]),
            reverse=True,
        )
        for next_node in successors:
            edge = (node, next_node)
            state_cap = state_caps.get(next_node, 1)
            edge_cap = edge_caps.get(edge, 1)
            if state_cap < 1:
                raise ValueError(f"Invalid state cap {state_cap} for {next_node}")
            if edge_cap < 1:
                raise ValueError(f"Invalid edge cap {edge_cap} for {edge}")
            if state_visits.get(next_node, 0) >= state_cap:
                continue
            if edge_visits.get(edge, 0) >= edge_cap:
                continue

            next_path = [*path, next_node]
            if len(next_path) > max_nodes:
                continue
            if next_node == end:
                walks.append(next_path)
                if len(walks) > max_walks:
                    raise RuntimeError(
                        "Assembler-bounded candidate population exceeded "
                        f"max_walks={max_walks}; enumeration was not silently truncated"
                    )
                continue

            distance = distance_to_end.get(next_node)
            if distance is None:
                continue
            remaining_slots = max_nodes - len(next_path)
            if distance > remaining_slots:
                continue

            next_state_visits = state_visits.copy()
            next_state_visits[next_node] += 1
            next_edge_visits = edge_visits.copy()
            next_edge_visits[edge] += 1
            stack.append(
                (next_node, next_path, next_state_visits, next_edge_visits)
            )

    walks.sort(key=format_oriented_path)
    return walks
