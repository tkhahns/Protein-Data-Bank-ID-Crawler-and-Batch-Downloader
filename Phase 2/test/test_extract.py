"""
This script contains unit tests for validating data extracted by extract methods 
against data retrieved using PDB's RESTful API. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "python -m unittest test/test_something".
To run all tests in the test directory, use the command "python -m unittest discover test" 
with the relevant flags (e.g. -v for higher verbosity).
"""

import unittest
import gemmi
from gemmi import cif, PolymerType
import os 
import glob
import requests
import requests_cache

from typing import NewType
import extract

MainData = NewType("MainData", tuple[str, str, str, str, str, str, str, str])
SubchainData = NewType("SubchainData", list) 

class TestExtractMethods(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.rootdir = "./database" # Location of .cif files
        self.POLYMER_TYPES = {"polypeptide(L)", "polypeptide(D)"}
        self.entry_ids = ['146D','178D','1A0A','1A0C','1A0F','1A0Q','1A1C','1A3I','1JKD','1LZH','1MM4','1PP3',
                          '1TGU','1XDF','2G9P','2HUM','3DSE','3IRL','3U7T','3UF8','4F5S','4W2P','5QB9','5SON',
                          '5U5C','5YII','6C6W','6FFL','6I06','6J4B','7A16','7H4H','7QTR','8E17','8FP7','8UHO','9B7F']
        requests_cache.install_cache('api_cache', expire_after=1800) # cache data for repeated calls 

    def read_file(self, entry_id):
        """Read the structure and document for a given entry id."""
        path = os.path.join(self.rootdir, f'*{entry_id.lower()}*')
        file_path = glob.glob(path)[0]

        struct = gemmi.read_structure(file_path)
        doc = cif.read(file_path)

        return struct, doc
    
    def make_request(self, url):
        try:
            response = requests.get(url, timeout=3)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print("Error occurred: ", e)

    def fetch_main_data(self, entry_id: str) -> MainData:
        """Fetch main table data for a given entry id (space_group, z_value, 
        a, b, c, alpha, beta, gamma). Any data not found in the JSON records are set to None."""

        try:
            main_data = self.make_request(f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}")

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

        except KeyError as e: 
            print(f"Error fetching data for entry {entry_id}, {e}")
        
    def fetch_subchain_data(self, entry_id: str) -> SubchainData: 
        """Fetch subchains sequence data from PDB for a given entry id. 
        Any data not found in the JSON records are set to None."""
        try:
            main_data = self.make_request(f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}") 

            polymer_entity_ids = main_data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
            data = []
            if polymer_entity_ids:
                for entity_id in polymer_entity_ids:
                    entity_data = self.make_request(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}")
                    polymer_type = entity_data.get('entity_poly', {}).get('type', None)
                    if polymer_type in self.POLYMER_TYPES: 
                        subchains = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('asym_ids', [])
                        for subchain in subchains:
                            seq = entity_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', None)
                            data.append(seq)
            return data
        except KeyError as e: 
            print(f"Error fetching data for entry {entry_id}, {e}")
        
        # subchain_ids: _rcsb_polymer_entity_container_identifiers.asym_ids
        # chain_ids: _rcsb_polymer_entity_container_identifiers.auth_asym_ids
    
    def test_insert_into_main_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                struct, doc = self.read_file(entry_id)
                result = extract.insert_into_main_table(struct, doc)[0][4:]
                expected = self.fetch_main_data(entry_id)

                self.assertEqual(result, expected)

    def test_subchain_sequence(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                struct, doc = self.read_file(entry_id)     
                result = extract.insert_into_subchain_table(struct, doc)
                result_seqs = [subchain[4] for subchain in result]
                #print(f'result seqs: {result_seqs}')
                expected_seqs = self.fetch_subchain_data(entry_id)
                #print(f'expected seqs: {expected_seqs}')
                self.assertCountEqual(result_seqs, expected_seqs)

    def tearDown(self):
        pass
    
if __name__ == "__main__":
    unittest.main()