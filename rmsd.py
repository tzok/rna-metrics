import sys

import numpy as np
from Bio import PDB


def calculate_rmsd(structure1_str, structure2_str):
    parser = PDB.PDBParser(QUIET=True)
    structure1 = parser.get_structure("structure1", structure1_str)
    structure2 = parser.get_structure("structure2", structure2_str)
    atoms1 = []
    atoms2 = []

    for model in structure1:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "P":  # Tylko atomy fosforu
                        atoms1.append(atom)

    for model in structure2:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "P":  # Tylko atomy fosforu
                        atoms2.append(atom)

    if len(atoms1) != len(atoms2):
        raise ValueError("Phosphorus atoms count mismatch")

    coords1 = np.array([atom.get_coord() for atom in atoms1])
    coords2 = np.array([atom.get_coord() for atom in atoms2])

    diff = coords1 - coords2
    rmsd = np.sqrt(np.sum(diff**2) / len(atoms1))
    return rmsd


def main(pdb_file1, pdb_file2):
    parser = PDB.PDBParser(QUIET=True)

    try:
        structure1 = parser.get_structure("structure1", pdb_file1)
        structure2 = parser.get_structure("structure2", pdb_file2)

        rmsd = calculate_rmsd(structure1, structure2)
        print(f"RMSD: {rmsd:.4f} Å")
    except Exception as e:
        raise RuntimeError(f"Nie można obliczyć RMSD: {str(e)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Użycie: python script.py <plik_pdb1> <plik_pdb2>")
        sys.exit(1)

    pdb_file1 = sys.argv[1]
    pdb_file2 = sys.argv[2]

    main(pdb_file1, pdb_file2)
