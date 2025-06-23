#! /usr/bin/env python
import sys
import tempfile
import os
import shutil
import subprocess

import numpy as np
from Bio.PDB import PDBParser
from rnapolis.unifier import main as unifier_main


def unify_structures(reference_pdb, model_pdb):
    """
    Unify two PDB structures using the rnapolis unifier.
    
    :param reference_pdb: Path to reference PDB file
    :param model_pdb: Path to model PDB file
    :return: Tuple of (unified_ref_path, unified_model_path, temp_dir)
    """
    temp_dir = tempfile.mkdtemp()
    
    # Save original sys.argv and replace it for unifier
    original_argv = sys.argv
    sys.argv = [
        "unifier",
        "--output",
        temp_dir,
        "--format",
        "keep",
        reference_pdb,
        model_pdb,
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
    ref_base = os.path.splitext(os.path.basename(reference_pdb))[0]
    model_base = os.path.splitext(os.path.basename(model_pdb))[0]
    ref_ext = os.path.splitext(reference_pdb)[1]
    model_ext = os.path.splitext(model_pdb)[1]

    unified_ref = os.path.join(temp_dir, f"{ref_base}{ref_ext}")
    unified_model = os.path.join(temp_dir, f"{model_base}{model_ext}")
    
    return unified_ref, unified_model, temp_dir


def calculate_lddt(reference_structure, model_structure):
    """
    Calculate the local Distance Difference Test (lDDT) score.

    :param reference_structure: PDB structure of the reference
    :param model_structure: PDB structure of the model
    :return: lDDT score (0-1, where 1 is perfect agreement)
    """
    ref_atoms = list(reference_structure.get_atoms())
    model_atoms = list(model_structure.get_atoms())

    if len(ref_atoms) != len(model_atoms):
        raise ValueError(
            "Number of atoms in reference and model structures do not match"
        )

    ref_coords = np.array([atom.coord for atom in ref_atoms])
    model_coords = np.array([atom.coord for atom in model_atoms])

    # Calculate pairwise distances between all atoms efficiently
    # This creates a 2D array where entry (i,j) is the distance between atom i and atom j
    ref_distances = np.linalg.norm(ref_coords[:, None] - ref_coords, axis=2)
    model_distances = np.linalg.norm(model_coords[:, None] - model_coords, axis=2)

    # Create a mask for interactions: atoms within 5A and not on the diagonal (self-interaction)
    interaction_mask = (ref_distances <= 5) & ~np.eye(len(ref_atoms), dtype=bool)

    # Exclude interactions between atoms in the same residue
    for i, atom in enumerate(ref_atoms):
        # Create a boolean array: True if atoms are in different residues, False otherwise
        different_residue = ~np.array(
            [atom.parent.id == ref_atoms[j].parent.id for j in range(len(ref_atoms))]
        )
        # Update the interaction mask to only include atoms in different residues
        interaction_mask[i] &= different_residue

    thresholds = [0.5, 1, 2, 4]
    scores = []

    for threshold in thresholds:
        preserved_interactions = np.abs(ref_distances - model_distances) < threshold
        score = np.sum(preserved_interactions[interaction_mask]) / np.sum(
            interaction_mask
        )
        scores.append(score)

    return np.mean(scores)


def main(reference_pdb, model_pdb):
    # Unify structures
    unified_ref, unified_model, temp_dir = unify_structures(reference_pdb, model_pdb)
    
    try:
        # Parse the unified structures
        parser = PDBParser()
        reference_structure = parser.get_structure("reference", unified_ref)
        model_structure = parser.get_structure("model", unified_model)

        lddt_score = calculate_lddt(reference_structure, model_structure)
        print(f"{lddt_score:.4f}")
    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python lddt.py <reference_pdb> <model_pdb>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
