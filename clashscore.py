#! /usr/bin/env python
import sys
import time
from pathlib import Path

import requests
from bs4 import BeautifulSoup


def calculate_clashscore(pdb_file):
    """Calculate clashscore using MolProbity web service."""
    breakpoint()
    base_url = "http://molprobity.biochem.duke.edu"
    session = requests.Session()

    # Step 1: Get the initial page and extract MolProbSID
    response = session.get(f"{base_url}/")
    soup = BeautifulSoup(response.text, "html.parser")
    molprobsid = soup.find("input", {"name": "MolProbSID"})["value"]

    # Step 2: Upload the file
    response = session.get(f"{base_url}/index.php?MolProbSID={molprobsid}")
    soup = BeautifulSoup(response.text, "html.parser")
    upload_event_id = soup.find("input", {"name": "eventID"})["value"]

    with open(pdb_file, "rb") as f:
        upload_data = {
            "MolProbSID": molprobsid,
            "cmd": "Upload >",
            "eventID": upload_event_id,
            "fetchType": "pdb",
            "pdbCode": "",
            "uploadType": "pdb",
        }
        files = {"uploadFile": (Path(pdb_file).name, f)}
        response = session.post(f"{base_url}/index.php", data=upload_data, files=files)

    # Step 3: Wait for processing and follow meta refreshes
    while True:
        response = session.get(f"{base_url}/index.php?MolProbSID={molprobsid}")
        soup = BeautifulSoup(response.text, "html.parser")

        # Check for meta refresh
        meta_refresh = soup.find("meta", {"http-equiv": "refresh"})
        if meta_refresh:
            content = meta_refresh["content"]
            if "; URL=" in content:
                redirect_url = content.split("; URL=")[1]
                response = session.get(redirect_url)
                continue

        # Check for Continue button
        continue_button = soup.find("input", {"type": "submit", "value": "Continue >"})
        if continue_button:
            form = continue_button.find_parent("form")
            if form:
                event_id = form.find("input", {"name": "eventID"})["value"]
                continue_data = {
                    "MolProbSID": molprobsid,
                    "cmd": "Continue >",
                    "eventID": event_id,
                }
    session.post(f"{base_url}/index.php", data=continue_data)

    # Step 5: Wait and get the analysis link
    while True:
        response = session.get(f"{base_url}/index.php?MolProbSID={molprobsid}")
        if response.status_code == 200:
            break
        time.sleep(1)

    # Find and follow "Analyze geometry without all-atom contacts" link
    soup = BeautifulSoup(response.text, "html.parser")
    analyze_link = soup.find("a", string="Analyze geometry without all-atom contacts")
    if not analyze_link:
        return None

    session.get(analyze_link["href"])

    # Step 6: Run the analysis
    file_name = Path(pdb_file).stem  # Remove .pdb extension
    analysis_data = {
        "MolProbSID": molprobsid,
        "chartAltloc": "1",
        "chartClashlist": "1",
        "chartNotJustOut": "1",
        "cmd": "Run programs to perform these analyses >",
        "doCharts": "1",
        "eventID": "61",
        "kinBaseP": "1",
        "kinGeom": "1",
        "kinSuite": "1",
        "modelID": file_name,
    }
    session.post(f"{base_url}/index.php", data=analysis_data)

    # Step 7: Wait for results and parse
    while True:
        response = session.get(f"{base_url}/index.php?MolProbSID={molprobsid}")
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, "html.parser")
            # Find the cell containing "Clashscore, all atoms:"
            clash_cell = soup.find("td", string="Clashscore, all atoms:")
            if clash_cell:
                # Get the next td element which contains the value
                value_cell = clash_cell.find_next("td")
                if value_cell:
                    try:
                        return float(value_cell.text.strip())
                    except ValueError:
                        return None
            break
        time.sleep(1)

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
