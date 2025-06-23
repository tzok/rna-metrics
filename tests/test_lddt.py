import pytest
import os
import shutil
from lddt import calculate_lddt, unify_structures
from Bio.PDB import PDBParser


class TestLDDT:
    def setup_method(self):
        """Set up test fixtures with paths to test PDB files."""
        self.pdb1 = "tests/1ehz.pdb"
        self.pdb2 = "tests/1evv.pdb"

        # Verify test files exist
        assert os.path.exists(self.pdb1), f"Test file {self.pdb1} not found"
        assert os.path.exists(self.pdb2), f"Test file {self.pdb2} not found"

    def test_lddt_calculation(self):
        """Test lDDT calculation between two structures using unified structures."""
        # Unify structures
        unified_ref, unified_model, temp_dir = unify_structures(self.pdb1, self.pdb2)

        try:
            # Parse the unified structures
            parser = PDBParser()
            reference_structure = parser.get_structure("reference", unified_ref)
            model_structure = parser.get_structure("model", unified_model)

            lddt_score = calculate_lddt(reference_structure, model_structure)
            assert lddt_score == pytest.approx(0.9907, abs=1e-4), (
                "lDDT score should be 0.9907"
            )
        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir)
