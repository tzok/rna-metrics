#! /usr/bin/env python
import sys

from rnapolis.annotator import extract_base_interactions
from rnapolis.common import LeontisWesthof
from rnapolis.parser import read_3d_structure


def calculate_inf(interactions1, interactions2):
    """Calculate INF score between two sets of interactions."""
    set1 = set(interactions1)
    set2 = set(interactions2)
    
    tp = len(set1 & set2)
    fn = len(set1 - set2)
    fp = len(set2 - set1)
    
    if tp == 0:
        return 0.0
    
    return (tp / (tp + fn) * tp / (tp + fp)) ** 0.5


def extract_interactions(interactions):
    """Extract different types of base interactions."""
    canonical_pairs = []
    non_canonical_pairs = []
    stacking_pairs = []

    # Process base pairs
    for pair in interactions.basePairs:
        if pair.lw == LeontisWesthof.cWW:
            seq = f"{pair.nt1.name}-{pair.nt2.name}"
            if seq in ["A-U", "U-A", "G-C", "C-G", "G-U", "U-G"]:
                canonical_pairs.append((pair.nt1, pair.nt2, None))
            else:
                non_canonical_pairs.append((pair.nt1, pair.nt2, pair.lw))
        else:
            non_canonical_pairs.append((pair.nt1, pair.nt2, pair.lw))

    # Process stacking interactions
    for stack in interactions.stackings:
        stacking_pairs.append((stack.nt1, stack.nt2, None))

    return canonical_pairs, non_canonical_pairs, stacking_pairs


def process_structure(pdb_file):
    """Process PDB file to extract different types of interactions."""
    with open(pdb_file) as f:
        structure = read_3d_structure(f)
        interactions = extract_base_interactions(structure)
        return extract_interactions(interactions)


def main(pdb_file1, pdb_file2):
    canonical1, non_canonical1, stacking1 = process_structure(pdb_file1)
    canonical2, non_canonical2, stacking2 = process_structure(pdb_file2)

    print(f"Structure {pdb_file1}:")
    print(f"  Canonical pairs: {len(canonical1)}")
    print(f"  Non-canonical pairs: {len(non_canonical1)}")
    print(f"  Stacking interactions: {len(stacking1)}")

    print(f"\nStructure {pdb_file2}:")
    print(f"  Canonical pairs: {len(canonical2)}")
    print(f"  Non-canonical pairs: {len(non_canonical2)}")
    print(f"  Stacking interactions: {len(stacking2)}")

    # Calculate INF scores for each interaction type
    canonical_inf = calculate_inf(canonical1, canonical2)
    non_canonical_inf = calculate_inf(non_canonical1, non_canonical2)
    stacking_inf = calculate_inf(stacking1, stacking2)

    print("\nINF Scores:")
    print(f"  Canonical pairs: {canonical_inf:.3f}")
    print(f"  Non-canonical pairs: {non_canonical_inf:.3f}")
    print(f"  Stacking interactions: {stacking_inf:.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python inf.py <reference_pdb> <model_pdb>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
