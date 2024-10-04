#! /usr/bin/env python
import sys

import numpy as np
from Bio import PDB


def extract_phosphorus_atoms(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "P":
                        atoms.append(atom)
    return atoms


def calculate_rmsd(structure1_str, structure2_str):
    parser = PDB.PDBParser(QUIET=True)
    structure1 = parser.get_structure("structure1", structure1_str)
    structure2 = parser.get_structure("structure2", structure2_str)

    atoms1 = extract_phosphorus_atoms(structure1)
    atoms2 = extract_phosphorus_atoms(structure2)

    if len(atoms1) != len(atoms2):
        raise ValueError("Phosphorus atoms count mismatch")

    coords1 = np.array([atom.get_coord() for atom in atoms1])
    coords2 = np.array([atom.get_coord() for atom in atoms2])

    diff = coords1 - coords2
    rmsd = np.sqrt(np.sum(diff**2) / len(atoms1))
    return rmsd


def main(pdb_file1, pdb_file2):
    try:
        with open(pdb_file1) as f:
            with open(pdb_file2) as g:
                print(calculate_rmsd(f.read(), g.read()))
    except Exception as e:
        raise RuntimeError(f"Nie można obliczyć RMSD: {str(e)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsd.py <pdb1> <pdb2>")
        sys.exit(1)

    pdb_file1 = sys.argv[1]
    pdb_file2 = sys.argv[2]

    main(pdb_file1, pdb_file2)
