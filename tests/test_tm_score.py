import pytest
import os
from tm_score import calculate_tm_score


class TestTMScore:
    def setup_method(self):
        """Set up test fixtures with paths to test PDB files."""
        self.pdb1 = "tests/1ehz.pdb"
        self.pdb2 = "tests/1evv.pdb"
        
        # Verify test files exist
        assert os.path.exists(self.pdb1), f"Test file {self.pdb1} not found"
        assert os.path.exists(self.pdb2), f"Test file {self.pdb2} not found"

    def test_tm_score_calculation(self):
        """Test TM-Score calculation between two structures."""
        tm_score = calculate_tm_score(self.pdb1, self.pdb2)
        assert tm_score == pytest.approx(0.9605, abs=1e-4), "TM-Score should be 0.9605"
