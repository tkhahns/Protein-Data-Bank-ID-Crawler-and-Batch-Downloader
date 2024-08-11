"""
This script contains integration tests for validating data in the database against data extracted by extract methods. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "python -m unittest test/test_something"
To run all tests in the test directory, use the command "python -m unittest discover test" 
with the relevant flags (e.g. -v for higher verbosity)
"""

import unittest
import sqlite3
import os
import re
import glob
from tqdm import tqdm
import gemmi 
from gemmi import cif 

import commands
import extract 

class TestPdbMethods(unittest.TestCase):

    def setUp(self):
        self.rootdir = "./database" # Location of .cif files
        self.entry_ids = ['146D','178D','1A0A','1A0C','1A0F','1A0Q','1A1C','1A3I','1JKD','1LZH','1MM4','1PP3',
                          '1TGU','1XDF','2G9P','2HUM','3DSE','3IRL','3U7T','3UF8','4F5S','4W2P','5QB9','5SON',
                          '5U5C','5YII','6C6W','6FFL','6I06','6J4B','7A16','7H4H','7QTR','8E17','8FP7','8UHO','9B7F']
        self.con = sqlite3.connect(':memory:')
        self.cur = self.con.cursor()
        commands.init_database(self.cur)

        for subdir, dirs, files in tqdm(os.walk(self.rootdir)):
            for file in files:
                path = os.path.join(subdir, file)
                if re.search('.*\.cif.*', path):
                    commands.insert_file(self.cur, path, verbose=False)
            self.con.commit()
        
    def read_file(self, entry_id):
        """Read the structure and document for a given entry id."""
        path = os.path.join(self.rootdir, f'*{entry_id.lower()}*')
        file_path = glob.glob(path)[0]

        struct = gemmi.read_structure(file_path)
        doc = cif.read(file_path)

        return struct, doc
    
    def query_data(self, entry_id, table, col="*"):
        """Query data from the database for a given entry id, table and columns (default all)."""
        query = f"SELECT {col} FROM {table} WHERE entry_id = ?"
        res = self.cur.execute(query, (entry_id,))

        return res
    
    def validate_data(self, entry_id, table, extract_func):
        """Validate data in the database against data extracted using the given function."""
        struct, doc = self.read_file(entry_id)
        result = self.query_data(entry_id, table).fetchall()
        expected = extract_func(struct, doc)

        self.assertCountEqual(result, expected)

    def test_main_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "main", extract.insert_into_main_table)
    
    def test_entity_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "entities", extract.insert_into_entity_table)

    def test_subchain_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "subchains", extract.insert_into_subchain_table)
                            
    def test_helix_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "helices", extract.insert_into_helix_table)

    def test_sheet_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "sheets", extract.insert_into_sheet_table)

    def test_strand_table(self):
        for entry_id in self.entry_ids:
            with self.subTest(entry_id=entry_id):
                self.validate_data(entry_id, "strands", extract.insert_into_strand_table)

    def tearDown(self):
        self.con.close()

if __name__ == "__main__":
    unittest.main()