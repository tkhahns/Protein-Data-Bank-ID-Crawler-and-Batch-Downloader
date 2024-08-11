"""
This script contains integration tests for validating data extracted by extract methods 
against data fetched from PDB's RESTful API. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
import gemmi
from gemmi import cif
import os
import glob
import requests
import requests_cache
import extract
from typing import NewType

MainData = NewType("MainData", tuple[str, str, str, str, str, str, str, str])
SubchainData = NewType("SubchainData", list[tuple[str, str, str, str, str, int]])

entry_ids = [
    '146D','178D','1A0A','1A0C','1A0F','1A0Q','1A1C','1A3I','1JKD','1LZH','1MM4','1PP3',
    '1TGU','1XDF','2G9P','2HUM','3DSE','3IRL','3U7T','3UF8','4F5S','4W2P','5QB9','5SON',
    '5U5C','5YII','6C6W','6FFL','6I06','6J4B','7A16','7H4H','7QTR','8E17','8FP7','8UHO','9B7F'
]

# Cache data for repeated calls
requests_cache.install_cache('api_cache', expire_after=1800)  

@pytest.fixture(scope="module")
def rootdir() -> str:
    """Fixture to provide the root directory of .cif files."""
    return "./database"  

@pytest.fixture(scope="module")
def polymer_types() -> set[str]:
    return {"polypeptide(L)", "polypeptide(D)"}

def read_file(rootdir: str, entry_id: str) -> tuple[gemmi.Structure, cif.Document]:
    """Read the structure and document for a given entry id."""
    path = os.path.join(rootdir, f'*{entry_id.lower()}*')
    file_path = glob.glob(path)[0]

    struct = gemmi.read_structure(file_path)
    doc = cif.read(file_path)

    return struct, doc

def make_request(obj: str, entry_id: str = "", entity_id: str = ""):
    """Make a request to retrieve PDB data through an endpoint, which differs according to the given data object (obj)."""
    try:
        urls = {
            "entry": f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
            "branched_entity": f"https://data.rcsb.org/rest/v1/core/branched_entity/{entry_id}/{entity_id}",
            "polymer_entity": f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}",
            "non_polymer_entity": f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{entry_id}/{entity_id}"
        }
        response = requests.get(urls[obj], timeout=3)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Request failed for {entry_id}: ", e)
        return {}

def fetch_main_data(entry_id: str) -> MainData:
    """Fetch main table data for a given entry id (space_group, z_value, a, b, c, alpha, beta, gamma). 
    Any data not found in the JSON records are set to None."""

    main_data = make_request("entry", entry_id=entry_id)

    spacegroup = main_data.get('symmetry', {}).get('space_group_name_hm', None)
    cell_data = main_data.get('cell', {})
    z_value = cell_data.get('zpdb', None)
    a = cell_data.get('length_a', None)
    b = cell_data.get('length_b', None)
    c = cell_data.get('length_c', None)
    alpha = cell_data.get('angle_alpha', None)
    beta = cell_data.get('angle_beta', None)
    gamma = cell_data.get('angle_gamma', None)

    return (spacegroup, z_value, a, b, c, alpha, beta, gamma)

def fetch_subchain_data(entry_id: str, polymer_types: set[str]) -> SubchainData:
    """Fetch subchains sequence data for a given entry id. Any data not found in the JSON records are set to None."""
    main_data = make_request("entry", entry_id=entry_id)
    polymer_entity_ids = main_data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
    data = []
    if polymer_entity_ids:
        for entity_id in polymer_entity_ids:
            entity_data = make_request("polymer_entity", entry_id=entry_id, entity_id=entity_id)
            polymer_type = entity_data.get('entity_poly', {}).get('type', None)
            if polymer_type in polymer_types:
                subchains = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('asym_ids', [])
                for index, subchain in enumerate(subchains):
                    chains = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', [])
                    chain_id = chains[index] # auth_asym_ids should be present in all right??
                    seq = entity_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', None)
                    length = entity_data.get('entity_poly', {}).get('rcsb_sample_sequence_length', None)
                    data.append((entry_id, entity_id, subchain, chain_id, seq, length))
    return data

@pytest.mark.parametrize("entry_id", entry_ids)
def test_insert_into_main_table(entry_id, rootdir):
    struct, doc = read_file(rootdir, entry_id)
    result = extract.insert_into_main_table(struct, doc)[0][4:]
    expected = fetch_main_data(entry_id)
        
    assert result == expected

@pytest.mark.parametrize("entry_id", entry_ids)
def test_insert_into_subchain_table(entry_id, rootdir, polymer_types):
    # Skip test if the entry does not contain protein entities, hence no subchains to test
    if entry_id in ['146D', '178D', '3IRL']: 
        pytest.skip(f"Entry {entry_id} does not contain protein entities")  
    struct, doc = read_file(rootdir, entry_id)
    result = extract.insert_into_subchain_table(struct, doc)
    result = [(subchain[0], subchain[1], subchain[2], subchain[3], subchain[4], subchain[7]) for subchain in result]  # all except start_pos, end_pos
    expected = fetch_subchain_data(entry_id, polymer_types)
    
    assert sorted(result) == sorted(expected)