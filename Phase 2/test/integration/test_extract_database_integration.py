"""
This script contains integration tests for validating data in the database against data extracted by extract methods. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/integration/test_something.py"
To run all tests in the test directory, use the command "pytest test/" 
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
import sqlite3
import os
import re
import glob
from tqdm import tqdm
import gemmi
from gemmi import cif
from typing import Generator, Callable, TypeVarTuple

import commands
import extract
from polymer_sequence import PolymerSequence

entry_ids = [
    '146D','178D','1A0A','1A0C','1A0F','1A0Q','1A1C','1A3I','1JKD','1LZH','1MM4','1PP3',
    '1TGU','1XDF','2G9P','2HUM','3DSE','3IRL','3U7T','3UF8','4F5S','4W2P','5QB9','5SON',
    '5U5C','5YII','6C6W','6FFL','6I06','6J4B','7A16','7H4H','7QTR','8E17','8FP7','8UHO','9B7F'
]

AttributeTypes = TypeVarTuple('AttributeTypes')

@pytest.fixture(scope="module")
def rootdir() -> str:
    return "./database"  # Location of .cif files

@pytest.fixture(scope="module")
def setup_database(rootdir: str) -> Generator[sqlite3.Cursor, None, None]:
    """Fixture to set up an in-memory database and insert PDB data for the given entry ids."""
    # Set up database
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    commands.init_database(cur)

    for subdir, dirs, files in tqdm(os.walk(rootdir)):
        for file in files:
            path = os.path.join(subdir, file)
            if re.search('.*\.cif.*', path):
                commands.insert_file(cur, path, verbose=False)
        con.commit()

    yield cur
    # Close database connection after all tests are run
    con.close() 

def read_file(rootdir: str, entry_id: str) -> tuple[gemmi.Structure, cif.Document, PolymerSequence]:
    """Read the structure, document and polymer sequence for a given entry id."""
    path = os.path.join(rootdir, f'*{entry_id.lower()}*')
    file_path = glob.glob(path)[0]

    struct = gemmi.read_structure(file_path)
    doc = cif.read(file_path)
    sequence = PolymerSequence(doc)

    return struct, doc, sequence

def query_data(cur: sqlite3.Cursor, entry_id: str, table: str, col: str = "*") -> sqlite3.Cursor:
    """Query data from the database for a given entry id, table and columns (default all)."""
    query = f"SELECT {col} FROM {table} WHERE entry_id = ?"
    res = cur.execute(query, (entry_id,))
    return res

def validate_data(entry_id: str, rootdir: str, cur: sqlite3.Cursor, table: str, extract_func: Callable[[gemmi.Structure, cif.Document], list[tuple[*AttributeTypes]]]):
    """Validate data in the database against data extracted using the given function."""
    struct, doc, sequence = read_file(rootdir, entry_id)
    result = query_data(cur, entry_id, table).fetchall()
    expected = extract_func(struct, doc, sequence)
    assert sorted(result) == sorted(expected)

@pytest.mark.parametrize("entry_id", entry_ids)  # run test for each entry_id
def test_main_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "main", extract.insert_into_main_table)

@pytest.mark.xfail(reason="extracted data are strings, not floats")
@pytest.mark.parametrize("entry_id", entry_ids)  
def test_experimental_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "experimental", extract.insert_into_experimental_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_entity_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "entities", extract.insert_into_entity_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_subchain_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "subchains", extract.insert_into_subchain_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_chain_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "chains", extract.insert_into_chain_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_helix_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "helices", extract.insert_into_helix_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_sheet_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "sheets", extract.insert_into_sheet_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_strand_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "strands", extract.insert_into_strand_table)

@pytest.mark.parametrize("entry_id", entry_ids)  
def test_coil_table(entry_id: str, rootdir: str, setup_database: sqlite3.Cursor):
    cur = setup_database
    validate_data(entry_id, rootdir, cur, "coils", extract.insert_into_coil_table)