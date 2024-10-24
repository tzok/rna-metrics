#! /usr/bin/env python
import sys

from rnapolis.annotator import extract_base_interactions
from rnapolis.parser import read_3d_structure


def extract_canonical_pairs(interactions):
    """Extract canonical base pairs (cWW A-U, G-C, G-U) from interactions."""
    canonical_pairs = []
    for pair in interactions.basePairs:
        if pair.leontisWesthof == "cWW":
            seq = f"{pair.nt1.residueName}-{pair.nt2.residueName}"
            if seq in ["A-U", "U-A", "G-C", "C-G", "G-U", "U-G"]:
                canonical_pairs.append(pair)
    return canonical_pairs


def main(pdb_file1, pdb_file2):
    structure1 = read_3d_structure(pdb_file1)
    structure2 = read_3d_structure(pdb_file2)
    interactions1 = extract_base_interactions(structure1)
    interactions2 = extract_base_interactions(structure2)
    
    canonical_pairs1 = extract_canonical_pairs(interactions1)
    canonical_pairs2 = extract_canonical_pairs(interactions2)
    
    print(f"Number of canonical pairs in {pdb_file1}: {len(canonical_pairs1)}")
    print(f"Number of canonical pairs in {pdb_file2}: {len(canonical_pairs2)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python inf.py <reference_pdb> <model_pdb>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
