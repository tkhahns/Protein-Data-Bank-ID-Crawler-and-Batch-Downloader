"""
This script contains unit tests for testing methods in table.py.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""
import pytest 
from unittest.mock import MagicMock
from table import Table

def test_table_initialisation():
    mock_attributes = MagicMock()
    mock_extractor = MagicMock(return_value=[("test_id", "data1")])
    test_table = Table("test_table", mock_attributes, mock_extractor)

    assert test_table.name == "test_table"
    assert test_table.attributes == mock_attributes
    assert test_table.extractor == mock_extractor

def test_table_initialisation_missing_arguments():
    with pytest.raises(TypeError):
        Table()

@pytest.mark.xfail(reason="attributes is not an iterable")
def test_attributes_string(test_table):
    test_table.attributes_string()

def test_create_table(test_table):
    expected = "CREATE TABLE IF NOT EXISTS test_table (id VARCHAR, a FLOAT,\
        PRIMARY KEY(id, a), FOREIGN KEY (id) REFERENCES\
                main(id))"
    result = test_table.create_table()

    assert expected == result

def test_retrieve_default_columns(test_table):
    expected = "SELECT * FROM test_table"
    result = test_table.retrieve()

    assert expected == result

def test_retreive_columns_specified(test_table):
    expected = "SELECT col1, col2 FROM test_table"
    test_columns = ("col1", "col2")
    result = test_table.retrieve(test_columns)

    assert expected == result

def test_retrieve_invalid_columns(test_table):
    test_columns = (1, 2)  # invalid as it is a tuple of ints
    with pytest.raises(TypeError):
        test_table.retrieve(test_columns)

def test_extract_data(test_table):
    mock_struct = MagicMock()
    mock_doc = MagicMock()
    mock_polymer_sequence = MagicMock()

    expected = [("test_id", "data1")]
    result = test_table.extract_data(mock_struct, mock_doc, mock_polymer_sequence)

    # check extractor function called with correct arguments 
    test_table.extractor.assert_called_with(mock_struct, mock_doc, mock_polymer_sequence)
    assert expected == result

def test_insert_row(test_table):
    test_data = ("test_id", "data1")
    expected = "INSERT INTO test_table VALUES(?, ?)"
    result = test_table.insert_row(test_data)

    assert expected == result

def test_insert_row_invalid_data(test_table):
    invalid_test_data = 1
    with pytest.raises(TypeError):
        test_table.insert_row(invalid_test_data)

@pytest.mark.xfail(reason="output differs from expected")
def test_update_row(test_table):
    test_data = {"col1": "value1", "col2": "value2"}
    test_primary_key_values = ["1"]
    expected = "UPDATE test_table SET col1 = value1, col2 = value2 WHERE entry_id = 1"
    result = test_table.update_row(test_data, test_primary_key_values)

    assert expected == result, f"Expected: {expected}, Got: {result}"
    # check that both methods were called once
    test_table.attributes.match_columns.assert_called_once_with(test_data)
    test_table.attributes.match_primary_keys.assert_called_once_with(test_primary_key_values)

def test_update_row_invalid_data(test_table):
    """
    Test that a ValueError is raised if match_primary_keys method raises a ValueError. 
    """
    test_table.attributes.match_columns.side_effect = ValueError("Number of values given does not match number of columns")
    test_data = {"col1": "value1", "col2": "value2"}
    test_primary_key_values = ["1"]
    
    with pytest.raises(ValueError, match="Number of values given does not match number of columns"):
        test_table.update_row(test_data, test_primary_key_values)
    
    test_table.attributes.match_columns.assert_called_once_with(test_data)
    # check that match_primary_keys was not called 
    test_table.attributes.match_primary_keys.assert_not_called()

def test_update_row_invalid_primary_keys(test_table):
    """
    Test that a ValueError is raised if match_primary_keys method raises a ValueError. 
    """
    test_table.attributes.match_primary_keys.side_effect = ValueError("Number of values given does not match number of keys")
    test_data = {"col1": "value1", "col2": "value2"}
    test_primary_key_values = ["1"]

    with pytest.raises(ValueError, match="Number of values given does not match number of keys"):
        test_table.update_row(test_data, test_primary_key_values)

    # check that both methods were called once
    test_table.attributes.match_columns.assert_called_once_with(test_data)
    test_table.attributes.match_primary_keys.assert_called_once_with(test_primary_key_values)
