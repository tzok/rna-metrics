import pytest
import os
import sys
import tempfile
from lddt import calculate_lddt, main as lddt_main
from Bio.PDB import PDBParser
from rnapolis.unifier import main as unifier_main
from io import StringIO


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
        # Create temporary directory for unified structures
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save original sys.argv and replace it for unifier
            original_argv = sys.argv
            sys.argv = [
                "unifier",
                "--output",
                temp_dir,
                "--format",
                "keep",
                self.pdb1,
                self.pdb2,
            ]

            try:
                unifier_main()
            except SystemExit:
                # unifier_main() calls sys.exit(), catch it to continue
                pass
            finally:
                # Restore original sys.argv
                sys.argv = original_argv

            # Get the unified file paths
            ref_base = os.path.splitext(os.path.basename(self.pdb1))[0]
            model_base = os.path.splitext(os.path.basename(self.pdb2))[0]
            ref_ext = os.path.splitext(self.pdb1)[1]
            model_ext = os.path.splitext(self.pdb2)[1]

            unified_ref = os.path.join(temp_dir, f"{ref_base}{ref_ext}")
            unified_model = os.path.join(temp_dir, f"{model_base}{model_ext}")

            # Parse the unified structures
            parser = PDBParser()
            reference_structure = parser.get_structure("reference", unified_ref)
            model_structure = parser.get_structure("model", unified_model)

            lddt_score = calculate_lddt(reference_structure, model_structure)
            assert lddt_score == pytest.approx(0.9907, abs=1e-4), (
                "lDDT score should be 0.9907"
            )
