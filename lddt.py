import numpy as np
from Bio.PDB import PDBParser
from rmsd import extract_phosphorus_atoms


def calculate_lddt(reference_structure, model_structure, cutoffs=[0.5, 1, 2, 4]):
    """
    Calculate the local Distance Difference Test (lDDT) score.

    :param reference_structure: PDB structure of the reference
    :param model_structure: PDB structure of the model
    :param cutoffs: Distance cutoffs for lDDT calculation (default: [0.5, 1, 2, 4])
    :return: lDDT score (0-1, where 1 is perfect agreement)
    """
    ref_atoms = extract_phosphorus_atoms(reference_structure)
    model_atoms = extract_phosphorus_atoms(model_structure)

    if len(ref_atoms) != len(model_atoms):
        raise ValueError(
            "Number of atoms in reference and model structures do not match"
        )

    n_atoms = len(ref_atoms)
    ref_coords = np.array([atom.coord for atom in ref_atoms])
    model_coords = np.array([atom.coord for atom in model_atoms])

    ref_distances = np.linalg.norm(ref_coords[:, None] - ref_coords, axis=2)
    model_distances = np.linalg.norm(model_coords[:, None] - model_coords, axis=2)

    scores = []
    for cutoff in cutoffs:
        ref_contacts = ref_distances < cutoff
        preserved_contacts = np.abs(ref_distances - model_distances) < 0.5
        score = np.sum(ref_contacts & preserved_contacts) / np.sum(ref_contacts)
        scores.append(score)

    return np.mean(scores)


def main(reference_pdb, model_pdb):
    parser = PDBParser()
    reference_structure = parser.get_structure("reference", reference_pdb)
    model_structure = parser.get_structure("model", model_pdb)

    lddt_score = calculate_lddt(reference_structure, model_structure)
    print(f"lDDT score: {lddt_score:.4f}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python lddt.py <reference_pdb> <model_pdb>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
