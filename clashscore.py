#! /usr/bin/env python
import sys
import urllib3

# Suppress InsecureRequestWarning
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import time
from pathlib import Path

import requests


def calculate_clashscore(pdb_file):
    """Calculate clashscore using MolProbity web service."""
    # MolProbity upload URL
    upload_url = "https://molprobity.biochem.duke.edu/molprobity/upload"
    results_url = "https://molprobity.biochem.duke.edu/molprobity/cmd"

    # Read PDB file
    with open(pdb_file, "rb") as f:
        files = {"file": (Path(pdb_file).name, f)}

        # Upload the structure
        response = requests.post(upload_url, files=files, verify=False)
        if not response.ok:
            return None

        # Get the session ID from response
        session_id = response.cookies.get("JSESSIONID")
        if not session_id:
            return None

        # Prepare cookies for next request
        cookies = {"JSESSIONID": session_id}

        # Request clashscore calculation
        params = {"command": "clashscore", "input": Path(pdb_file).name}

        # Send request and wait for results
        max_attempts = 10
        for _ in range(max_attempts):
            response = requests.get(results_url, params=params, cookies=cookies, verify=False)
            if response.ok and "clashscore =" in response.text:
                # Extract clashscore from response
                for line in response.text.split("\n"):
                    if "clashscore =" in line:
                        try:
                            return float(line.split("=")[1].strip())
                        except (IndexError, ValueError):
                            return None
            time.sleep(1)  # Wait before retrying

        return None


def main(pdb_file):
    score = calculate_clashscore(pdb_file)
    if score is not None:
        print(f"{score:.4f}")
    else:
        print("nan")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clashscore.py <pdb_file>")
        sys.exit(1)
    main(sys.argv[1])
