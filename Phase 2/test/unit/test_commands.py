"""
This script contains unit tests for testing methods in commands.py.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""
import pytest
from unittest.mock import patch, call, MagicMock
import gemmi

import table
import commands 

TEST_FILE_PATH = "test_path/file.cif"
TEST_DATA = ('1A00', 'data1', 'data2')
TEST_STATEMENT = "INSERT INTO main VALUES(?, ?, ?)"

def test_init_database(mock_cursor):
    mock_table_1 = MagicMock(spec=table.Table)
    mock_statement_1 = "CREATE TABLE IF NOT EXISTS \
        test_table (id VARCHAR, PRIMARY KEY (id))"
    mock_table_1.create_table.return_value = mock_statement_1

    mock_table_2 = MagicMock(spec=table.Table)
    
    mock_statement_2 = "CREATE TABLE IF NOT EXISTS \
        test_table (id VARCHAR, PRIMARY KEY (id), \
            FOREIGN KEY (id) REFERENCES main (id))" 
    mock_table_2.create_table.return_value = mock_statement_2

    mock_table_schemas = [mock_table_1, mock_table_2]
    with patch('commands.table_schemas', mock_table_schemas):
        commands.init_database(mock_cursor)

        mock_cursor.execute.assert_any_call(mock_statement_1)
        mock_cursor.execute.assert_any_call(mock_statement_2)

        assert mock_cursor.execute.call_count == 2


def test_init_database_empty_table_schemas(mock_cursor):
    """
    Test that the execute method is not called when table_schemas
    is an empty list.
    """
    with patch('commands.table_schemas', []):
        commands.init_database(mock_cursor)

        mock_cursor.execute.assert_not_called()


def test_init_database_missing_create_table_method(mock_cursor, mock_table):
    """
    Test that an attribute error is raised if the create_table
    method is missing.
    """
    del mock_table.create_table  # delete create_table method

    mock_table_schemas = [mock_table]
    with patch('commands.table_schemas', mock_table_schemas):
        with pytest.raises(AttributeError):
            commands.init_database(mock_cursor)


@patch('gemmi.cif.read')
@patch('commands.PolymerSequence')
def test_insert_file_entry_not_in_main_table(mock_cif_read, mock_polymer_seq, mock_structure, mock_table, mock_cursor, capsys):
    """
    Test that data is inserted into the table when the entry
    is not found in the main table. 
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_cif_read.return_value = MagicMock()
        mock_polymer_seq.return_value = MagicMock()

        mock_cursor.execute.return_value.fetchone.return_value = None

        mock_table_schemas = [mock_table]
        with patch('commands.table_schemas', mock_table_schemas):
            commands.insert_file(mock_cursor, TEST_FILE_PATH)
            # check file_path was printed
            assert TEST_FILE_PATH in capsys.readouterr().out  
            expected_calls = [call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                              call.execute().fetchone(),
                              call.execute(TEST_STATEMENT, TEST_DATA)]
            mock_cursor.assert_has_calls(expected_calls)
        

@pytest.mark.xfail(reason="unable to access struct variable")
@patch('gemmi.read_structure')
def test_insert_file_exception_gemmi_read_failure(mock_gemmi_read, mock_cursor, capsys):
    """
    Test that the exception is caught when an error occurs when
    reading the gemmi Structure.
    """
    mock_gemmi_read.side_effect = Exception("Error reading structure")
    commands.insert_file(mock_cursor, TEST_FILE_PATH)
    captured = capsys.readouterr()

    # check file name, structure name and error are printed
    assert TEST_FILE_PATH in captured.out
    assert "mock_name" in captured.out 
    assert "Error reading structure" in captured.out


@patch('gemmi.cif.read')
@patch('commands.PolymerSequence')
def test_insert_file_exception_cif_read_failure(mock_cif_read, mock_polymer_seq, mock_structure, mock_cursor, capsys):
    """
    Test that the exception is caught when an error occurs when
    reading the cif Document.
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_cif_read.side_effect = Exception("Error reading document")
        mock_polymer_seq.return_value = MagicMock()
        commands.insert_file(mock_cursor, TEST_FILE_PATH)
        captured = capsys.readouterr()
        
        # check file name, structure name and error are printed
        assert TEST_FILE_PATH in captured.out
        assert "mock_name" in captured.out 
        assert "Error reading document" in captured.out 


@patch('gemmi.cif.read')
@patch('commands.PolymerSequence')
def test_insert_file_entry_exists_in_main_table(mock_cif_read, mock_polymer_seq, mock_structure, mock_table, mock_cursor):
    """
    Test that data is not inserted into the table when the entry 
    already exists in the main table.
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_cif_read.return_value = MagicMock()
        mock_polymer_seq.return_value = MagicMock()

        mock_cursor.execute.return_value.fetchone.return_value = ('1A00', )

        mock_table_schemas = [mock_table]
        with patch('commands.table_schemas', mock_table_schemas):
            commands.insert_file(mock_cursor, TEST_FILE_PATH)
            expected_calls = [call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                              call.execute().fetchone()]
            mock_cursor.assert_has_calls(expected_calls)
            # check cur.execute(statement, data) is not called
            assert call(TEST_STATEMENT, TEST_DATA) not in mock_cursor.execute.call_args_list
