"""
This script contains unit tests for testing methods in database.py.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, call, MagicMock, PropertyMock
import gemmi
from sqlite3 import OperationalError

import database

TEST_TABLE_NAME = "main"
TEST_ENTRY_ID = "1A00"
TEST_DATA = ('1A00', 'data1', 'data2')

def test_insert_into_table(mock_cursor): 
    expected_query = "INSERT INTO main VALUES(?, ?, ?)"
    database.insert_into_table(mock_cursor, TEST_TABLE_NAME, TEST_DATA)

    mock_cursor.execute.assert_called_with(expected_query, TEST_DATA)
    

def test_insert_into_table_empty_data(mock_cursor):
    empty_test_data = ()
    expected_query = "INSERT INTO main VALUES()"
    mock_cursor.execute.side_effect = OperationalError()
    
    with pytest.raises(OperationalError):
        database.insert_into_table(mock_cursor, TEST_TABLE_NAME, empty_test_data)
        mock_cursor.execute.assert_called_with(expected_query, empty_test_data)

def test_retrieve_from_table(mock_cursor):
    expected_query = "SELECT * FROM main WHERE entry_id = 1A00"
    mock_cursor.execute.return_value.fetchall.return_value = [TEST_DATA]

    result = database.retrieve_from_table(mock_cursor, TEST_TABLE_NAME, TEST_ENTRY_ID)
    
    assert result == [TEST_DATA]
    mock_cursor.execute.assert_called_once_with(expected_query)


def test_retrieve_from_table_entry_not_found(mock_cursor):
    expected_query = "SELECT * FROM main WHERE entry_id = 1A00"
    # empty list is returned since entry not found in table 
    mock_cursor.execute.return_value.fetchall.return_value = []

    result = database.retrieve_from_table(mock_cursor, TEST_TABLE_NAME, TEST_ENTRY_ID)
    
    assert result == []
    mock_cursor.execute.assert_called_once_with(expected_query)


