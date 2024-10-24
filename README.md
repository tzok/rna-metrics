# RNA Structure Analysis Tools

A collection of Python tools for analyzing RNA 3D structures and calculating various structural similarity metrics.

## Tools

### TM-score

Calculates the Template Modeling Score between two RNA structures using USalign.

Usage:

```bash
tm_score.py <reference_pdb> <model_pdb>
```

### RMSD (Root Mean Square Deviation)

Calculates RMSD between phosphorus atoms of two RNA structures after optimal superposition.

Usage:

```bash
rmsd.py <reference_pdb> <model_pdb>
```

### INF (Interaction Network Fidelity)

Measures the similarity of base-base interaction networks between two RNA structures.
Supports different interaction types: canonical base pairs, non-canonical base pairs, and stacking interactions.

Usage:

```bash
inf.py <reference_pdb> <model_pdb> [mode]
```

Mode options: canonical, non-canonical, stacking, all (default)

### MCQ (Mean of Circular Quantities)

Calculates the mean angular difference between torsion angles of two RNA structures.

Usage:

```bash
mcq.py <reference_pdb> <model_pdb>
```

### lDDT (local Distance Difference Test)

Evaluates local structure quality by comparing atomic distance patterns.

Usage:

```bash
lddt.py <reference_pdb> <model_pdb>
```

### Clashscore

Calculates atomic clashes using the MolProbity web service.

Usage:

```bash
clashscore.py <pdb_file>
```

### Torsion Angles

Calculates and compares RNA backbone torsion angles between two structures.

Usage:

```bash
torsion.py <pdb_file1> <pdb_file2>
```

## Docker Usage

Build the container:

```bash
docker build -t rna-tools .
```

Run a tool (example with TM-score):

```bash
docker run --rm -v $(pwd):/data rna-tools tm_score.py /data/reference.pdb /data/model.pdb
```

## Dependencies

- Python 3.12
- Biopython
- NumPy
- RNApolis
- BeautifulSoup4
- Requests

Install dependencies:

```bash
pip install -r requirements.txt
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
