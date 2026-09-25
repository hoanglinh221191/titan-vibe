#!/usr/bin/env python3
"""Focused deterministic tests for H9 graph descriptors."""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))

from analyze_h9_evidence_graph import build_graph_metrics  # noqa: E402


def expect_close(actual: float | None, expected: float, label: str) -> None:
    if actual is None or not math.isclose(actual, expected, rel_tol=1.0e-11, abs_tol=1.0e-12):
        raise AssertionError(f"{label}: expected {expected}, observed {actual}")


def brute_edge_connectivity(node_count: int, edges: list[tuple[int, int]]) -> int:
    best = len(edges)
    # Anchor node 0 on one side so complementary cuts are not duplicated.
    for mask in range(0, (1 << (node_count - 1)) - 1):
        crossing = 0
        for left, right in edges:
            left_side = left == 0 or bool(mask & (1 << (left - 1)))
            right_side = right == 0 or bool(mask & (1 << (right - 1)))
            crossing += left_side != right_side
        best = min(best, crossing)
    return best


def main() -> int:
    path3 = build_graph_metrics(3, [(0, 1), (1, 2)])
    assert path3.component_count == 1
    assert path3.edge_connectivity == 1
    assert path3.bridge_count == 2
    expect_close(path3.bridge_fraction, 1.0, "path3 bridge fraction")
    expect_close(path3.normalized_algebraic_connectivity, 1.0, "path3 lambda2")

    triangle = build_graph_metrics(3, [(0, 1), (1, 2), (0, 2)])
    assert triangle.edge_connectivity == 2
    assert triangle.bridge_count == 0
    expect_close(triangle.edge_density, 1.0, "triangle density")
    expect_close(
        triangle.normalized_algebraic_connectivity, 1.5, "triangle lambda2"
    )

    complete4 = build_graph_metrics(
        4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    )
    assert complete4.edge_connectivity == 3
    expect_close(
        complete4.normalized_algebraic_connectivity, 4.0 / 3.0, "K4 lambda2"
    )

    disconnected = build_graph_metrics(4, [(0, 1), (2, 3)])
    assert disconnected.component_count == 2
    assert disconnected.edge_connectivity == 0
    expect_close(
        disconnected.normalized_algebraic_connectivity, 0.0, "disconnected lambda2"
    )

    singleton = build_graph_metrics(1, [])
    assert singleton.edge_connectivity is None
    assert singleton.normalized_algebraic_connectivity is None
    assert singleton.bridge_fraction is None

    generator = random.Random(20260828)
    all_edges = [(left, right) for left in range(6) for right in range(left + 1, 6)]
    for _ in range(40):
        edges = [edge for edge in all_edges if generator.random() < 0.4]
        observed = build_graph_metrics(6, edges).edge_connectivity
        expected = brute_edge_connectivity(6, edges)
        assert observed == expected, (edges, observed, expected)

    print("PASS h9_evidence_graph_analysis_smoke graphs=5 random_edge_cuts=40")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
