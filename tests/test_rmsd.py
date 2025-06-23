import pytest
import os
from rmsd import calculate_rmsd


class TestRMSD:
    def setup_method(self):
        """Set up test fixtures with paths to test PDB files."""
        self.pdb1 = "tests/1ehz.pdb"
        self.pdb2 = "tests/1evv.pdb"
        
        # Verify test files exist
        assert os.path.exists(self.pdb1), f"Test file {self.pdb1} not found"
        assert os.path.exists(self.pdb2), f"Test file {self.pdb2} not found"

    def test_rmsd_calculation(self):
        """Test RMSD calculation between two structures."""
        with open(self.pdb1) as f1, open(self.pdb2) as f2:
            rmsd = calculate_rmsd(f1.read(), f2.read())
            assert rmsd == pytest.approx(0.5935, abs=1e-4), "RMSD score should be 0.5935"
