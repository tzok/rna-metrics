import pytest
import os
from mcq import calculate_mcq
from torsion import calculate_torsion_angles
from Bio import PDB


class TestMCQ:
    def setup_method(self):
        """Set up test fixtures with paths to test PDB files."""
        self.pdb1 = "tests/1ehz.pdb"
        self.pdb2 = "tests/1evv.pdb"

        # Verify test files exist
        assert os.path.exists(self.pdb1), f"Test file {self.pdb1} not found"
        assert os.path.exists(self.pdb2), f"Test file {self.pdb2} not found"

        # Load structures and calculate torsion angles
        parser = PDB.PDBParser(QUIET=True)
        structure1 = parser.get_structure("struct1", self.pdb1)
        structure2 = parser.get_structure("struct2", self.pdb2)

        self.angles1 = calculate_torsion_angles(structure1)
        self.angles2 = calculate_torsion_angles(structure2)

    def test_mcq_calculation(self):
        """Test MCQ calculation between two structures."""
        mcq = calculate_mcq(self.angles1, self.angles2)
        assert mcq == pytest.approx(9.5274, abs=1e-4), "MCQ score should be 9.5274"
