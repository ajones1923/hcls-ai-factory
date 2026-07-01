"""Tests for single-cell annotation logic (D1). Pure logic — no scanpy needed."""
from single_cell_compute import annotate_cluster, summarize_clusters


class TestAnnotate:
    def test_b_cells(self):
        r = annotate_cluster(["MS4A1", "CD79A", "ZZZ", "QQQ"])
        assert r["cell_type"] == "B cells" and r["marker_overlap"] == 2
        assert set(r["evidence"]) == {"MS4A1", "CD79A"}

    def test_monocytes(self):
        assert annotate_cluster(["CD14", "LYZ", "S100A8", "S100A9"])["cell_type"] == "CD14+ Monocytes"

    def test_unknown_when_no_marker_overlap(self):
        assert annotate_cluster(["FOO", "BAR", "BAZ"])["cell_type"] == "Unknown"


class TestSummarize:
    def test_summarize_annotates_each_cluster(self):
        clusters = [
            {"cluster": "0", "n_cells": 500, "top_genes": ["CD14", "LYZ", "S100A8"]},
            {"cluster": "1", "n_cells": 300, "top_genes": ["MS4A1", "CD79A"]},
            {"cluster": "2", "n_cells": 120, "top_genes": ["GNLY", "NKG7", "KLRD1"]},
        ]
        s = summarize_clusters(clusters)
        cts = {c["cluster"]: c["cell_type"] for c in s["clusters"]}
        assert cts == {"0": "CD14+ Monocytes", "1": "B cells", "2": "NK cells"}
        assert s["n_clusters"] == 3 and "NK cells" in s["cell_types"]
