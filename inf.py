#! /usr/bin/env python
import sys

from rnapolis.annotator import extract_base_interactions
from rnapolis.common import LeontisWesthof
from rnapolis.parser import read_3d_structure


def extract_canonical_pairs(interactions):
    """Extract canonical base pairs (cWW A-U, G-C, G-U) from interactions."""
    canonical_pairs = []
    for pair in interactions.basePairs:
        if pair.lw == LeontisWesthof.cWW:
            seq = f"{pair.nt1.name}-{pair.nt2.name}"
            if seq in ["A-U", "U-A", "G-C", "C-G", "G-U", "U-G"]:
                canonical_pairs.append(pair)
    return canonical_pairs


def process_structure(pdb_file):
    """Process PDB file to extract canonical base pairs."""
    with open(pdb_file) as f:
        structure = read_3d_structure(f)
        interactions = extract_base_interactions(structure)
        canonical_pairs = extract_canonical_pairs(interactions)
    return canonical_pairs


def main(pdb_file1, pdb_file2):
    canonical_pairs1 = process_structure(pdb_file1)
    canonical_pairs2 = process_structure(pdb_file2)

    print(f"Number of canonical pairs in {pdb_file1}: {len(canonical_pairs1)}")
    print(f"Number of canonical pairs in {pdb_file2}: {len(canonical_pairs2)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python inf.py <reference_pdb> <model_pdb>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
