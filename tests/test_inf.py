import pytest
import os
from inf import process_structure, calculate_inf


class TestINF:
    def setup_method(self):
        """Set up test fixtures with paths to test PDB files."""
        self.pdb1 = "tests/1ehz.pdb"
        self.pdb2 = "tests/1evv.pdb"

        # Verify test files exist
        assert os.path.exists(self.pdb1), f"Test file {self.pdb1} not found"
        assert os.path.exists(self.pdb2), f"Test file {self.pdb2} not found"

        # Process structures once for all tests
        self.canonical1, self.non_canonical1, self.stacking1 = process_structure(
            self.pdb1
        )
        self.canonical2, self.non_canonical2, self.stacking2 = process_structure(
            self.pdb2
        )

    def test_inf_canonical(self):
        """Test INF calculation for canonical base pairs."""
        score = calculate_inf(self.canonical1, self.canonical2)
        assert score == pytest.approx(1.0, abs=1e-4), (
            "Canonical INF score should be 1.0"
        )

    def test_inf_non_canonical(self):
        """Test INF calculation for non-canonical base pairs."""
        score = calculate_inf(self.non_canonical1, self.non_canonical2)
        assert score == pytest.approx(0.9375, abs=1e-4), (
            "Non-canonical INF score should be 0.9375"
        )

    def test_inf_stacking(self):
        """Test INF calculation for stacking interactions."""
        score = calculate_inf(self.stacking1, self.stacking2)
        assert score == pytest.approx(0.9506, abs=1e-4), (
            "Stacking INF score should be 0.9506"
        )

    def test_inf_all(self):
        """Test INF calculation for all interactions combined."""
        all_interactions1 = self.canonical1 + self.non_canonical1 + self.stacking1
        all_interactions2 = self.canonical2 + self.non_canonical2 + self.stacking2
        score = calculate_inf(all_interactions1, all_interactions2)
        assert score == pytest.approx(0.9570, abs=1e-4), (
            "All interactions INF score should be 0.9570"
        )
