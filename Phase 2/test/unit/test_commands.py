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


@patch("commands.insert_file")
@patch("gemmi.cif.read")
@patch("commands.PolymerSequence")
def test_check_file_entry_not_in_main_table(mock_polymer_seq, mock_cif_read, mock_insert_file, mock_structure, mock_cursor, capsys):
    """
    Test that data is inserted into the table when the entry is not found in the main table. 
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_doc = MagicMock()
        mock_cif_read.return_value = mock_doc
        mock_sequence = MagicMock()
        mock_polymer_seq.return_value = mock_sequence

        # no row in the main table with the entry ID 
        mock_cursor.execute.return_value.fetchone.return_value = None

        commands.check_file(mock_cursor, TEST_FILE_PATH)
        captured = capsys.readouterr().out 

        # check that file_path was printed
        assert "Checking " + TEST_FILE_PATH in captured 
        assert "Adding " + TEST_FILE_PATH in captured
        
        # check that cursor methods were called 
        expected_calls = [call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                            call.execute().fetchone()]
                            #call.execute(TEST_STATEMENT, TEST_DATA)]
        mock_cursor.assert_has_calls(expected_calls)  
        
        # check that insert_file was called 
        mock_insert_file.assert_called_once_with(mock_cursor, mock_structure, mock_doc, mock_sequence)
        

#@pytest.mark.xfail(reason="unable to access struct variable")
@patch("gemmi.read_structure")
def test_check_file_gemmi_read_failure(mock_gemmi_read, mock_cursor, capsys):
    """
    Test that the exception is caught when an error occurs when
    reading the gemmi Structure.
    """
    mock_gemmi_read.side_effect = Exception("Error reading structure")
    commands.check_file(mock_cursor, TEST_FILE_PATH)
    captured = capsys.readouterr()

    # check file name, structure name and error are printed
    assert "Checking " + TEST_FILE_PATH in captured.out
    assert "Error reading structure" in captured.out


@patch("gemmi.cif.read")
@patch("commands.PolymerSequence")
def test_check_file_cif_read_failure(mock_polymer_seq, mock_cif_read, mock_structure, mock_cursor, capsys):
    """
    Test that the exception is caught when an error occurs when
    reading the cif Document.
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_cif_read.side_effect = Exception("Error reading document")
        mock_polymer_seq.return_value = MagicMock()
        commands.check_file(mock_cursor, TEST_FILE_PATH)
        captured = capsys.readouterr()
        
        # check file name, structure name and error are printed
        assert "Checking " + TEST_FILE_PATH in captured.out
        assert "mock_name" in captured.out 
        assert "Error reading document" in captured.out 


@patch("commands.update_file")
@patch("gemmi.cif.read")
@patch("commands.PolymerSequence")
def test_check_file_entry_exists_needs_revision(mock_polymer_seq, mock_cif_read, mock_update_file, mock_structure, mock_table_schemas, mock_cursor, capsys):
    """
    Test that data is not inserted into the table when the entry 
    already exists in the main table and is revised when it is not up to date.
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_doc, mock_block = MagicMock(), MagicMock()
        mock_doc.sole_block.return_value = mock_block
        mock_block.find_value.return_value = "2000-12-31"  # mock revision_date
        mock_cif_read.return_value = mock_doc
        mock_sequence = MagicMock()
        mock_polymer_seq.return_value = mock_sequence

        # row in the main table with such entry ID exists, and is not up to date 
        mock_cursor.execute.return_value.fetchone.side_effect = [
            ('1A00', ),
            ("2000-12-01", )
        ]

        with patch('commands.table_schemas', mock_table_schemas):
            commands.check_file(mock_cursor, TEST_FILE_PATH)
            expected_calls = [
                call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                call.execute().fetchone(),
                call.execute("SELECT revision_date FROM main WHERE entry_id = '1A00'"),
                call.execute().fetchone()
            ]
            mock_cursor.assert_has_calls(expected_calls)
            captured = capsys.readouterr()
        
            # check that file name was printed
            assert "Checking " + TEST_FILE_PATH in captured.out
            assert "Updating " + TEST_FILE_PATH in captured.out
            # check that update file was called 
            mock_update_file.assert_called_once_with(mock_cursor, mock_structure, mock_doc, mock_sequence)


@patch("gemmi.cif.read")
@patch("commands.PolymerSequence")
def test_check_file_entry_exists_not_corrupted(mock_polymer_seq, mock_cif_read, mock_structure, mock_table_schemas, mock_cursor, capsys):
    """
    Test that data is not updated when the entry already exists in the main table, 
    is up to date and not corrupted. 
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_doc, mock_block = MagicMock(), MagicMock()
        mock_doc.sole_block.return_value = mock_block
        mock_block.find_value.return_value = "2000-12-31"  # mock revision_date
        mock_cif_read.return_value = mock_doc
        mock_sequence = MagicMock()
        mock_polymer_seq.return_value = mock_sequence

        # row in the main table with such entry ID exists, is up to date
        # and a row in the coils table with such entry ID exists 
        mock_cursor.execute.return_value.fetchone.side_effect = [
            ('1A00', ),
            ("2000-12-31", ),
            ('1A00', )
        ]

        with patch('commands.table_schemas', mock_table_schemas):
            commands.check_file(mock_cursor, TEST_FILE_PATH)
            expected_calls = [
                call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                call.execute().fetchone(),
                call.execute("SELECT revision_date FROM main WHERE entry_id = '1A00'"),
                call.execute().fetchone(), 
                call.execute("SELECT entry_id FROM coils WHERE entry_id = '1A00'"),
                call.execute().fetchone()
            ]
            mock_cursor.assert_has_calls(expected_calls)
            captured = capsys.readouterr()
        
            # check that file name was printed
            assert "Checking " + TEST_FILE_PATH in captured.out


@patch("commands.update_file")
@patch("gemmi.cif.read")
@patch("commands.PolymerSequence")
def test_check_file_entry_exists_is_corrupted(mock_polymer_seq, mock_cif_read, mock_update_file, mock_structure, mock_table_schemas, mock_cursor, capsys):
    """
    Test that data is updated when the entry already exists in the main table, 
    is up to date and corrupted. 
    """
    with patch.object(gemmi,'read_structure', return_value=mock_structure):
        mock_doc, mock_block = MagicMock(), MagicMock()
        mock_doc.sole_block.return_value = mock_block
        mock_block.find_value.return_value = "2000-12-31"  # mock revision_date
        mock_cif_read.return_value = mock_doc
        mock_sequence = MagicMock()
        mock_polymer_seq.return_value = mock_sequence

        # row in the main table with such entry ID exists, is up to date
        # but no row in the coils table with such entry ID exists 
        mock_cursor.execute.return_value.fetchone.side_effect = [
            ('1A00', ),
            ("2000-12-31", ),
            None
        ]

        with patch('commands.table_schemas', mock_table_schemas):
            commands.check_file(mock_cursor, TEST_FILE_PATH)
            expected_calls = [
                call.execute("SELECT entry_id FROM main WHERE entry_id = '1A00'"), 
                call.execute().fetchone(),
                call.execute("SELECT revision_date FROM main WHERE entry_id = '1A00'"),
                call.execute().fetchone(), 
                call.execute("SELECT entry_id FROM coils WHERE entry_id = '1A00'"),
                call.execute().fetchone()
            ]
            mock_cursor.assert_has_calls(expected_calls)
            captured = capsys.readouterr()
        
            # check that file name was printed
            assert "Checking " + TEST_FILE_PATH in captured.out
            assert "Data corrupted, fixing " + TEST_FILE_PATH in captured.out
            # check that update file was called 
            mock_update_file.assert_called_once_with(mock_cursor, mock_structure, mock_doc, mock_sequence)


def test_insert_file(mock_table_schemas, mock_cursor):
    with patch('commands.table_schemas', mock_table_schemas):
        commands.insert_file(mock_cursor, MagicMock(), MagicMock(), MagicMock())
        expected_calls = [
            call.execute("INSERT INTO main VALUES(?, ?, ?)", ('1A00', 'data1', 'data2')),
            call.execute("INSERT INTO coils VALUES(?, ?, ?)", ('1A00', 'data1', 'data2')),
            call.execute("INSERT INTO coils VALUES(?, ?, ?)", ('1A00', 'data3', 'data4'))
        ]
        
        mock_cursor.assert_has_calls(expected_calls)


def test_insert_file_no_data_to_extract(mock_table_schemas, mock_cursor):
    """
    Test valid insertion of data when no data is extracted for one of the tables. 
    """
    with patch('commands.table_schemas', mock_table_schemas):
        # extract_data returns an empty list for the coils table 
        mock_table_schemas[-1].extract_data.return_value = []

        commands.insert_file(mock_cursor, MagicMock(), MagicMock(), MagicMock())
        expected_calls = [
            call.execute("INSERT INTO main VALUES(?, ?, ?)", ('1A00', 'data1', 'data2'))
        ]
        
        mock_cursor.assert_has_calls(expected_calls)


def test_update_file(mock_table_schemas, mock_cursor, mock_structure):
    with patch('commands.table_schemas', mock_table_schemas):
        commands.update_file(mock_cursor, mock_structure, MagicMock(), MagicMock())
        expected_calls = [
            call.execute("DELETE FROM main WHERE entry_id = '1A00'"),
            call.execute("INSERT INTO main VALUES(?, ?, ?)", ('1A00', 'data1', 'data2')),
            call.execute("DELETE FROM coils WHERE entry_id = '1A00'"),
            call.execute("INSERT INTO coils VALUES(?, ?, ?)", ('1A00', 'data1', 'data2')),
            call.execute("INSERT INTO coils VALUES(?, ?, ?)", ('1A00', 'data3', 'data4'))
        ]
        
        mock_cursor.assert_has_calls(expected_calls)


def test_update_file_no_data_to_extract(mock_table_schemas, mock_cursor, mock_structure):
    """
    Test valid updating of data when no data is extracted for one of the tables. 
    """
    with patch('commands.table_schemas', mock_table_schemas):
         # extract_data returns an empty list for the coils table 
        mock_table_schemas[-1].extract_data.return_value = []

        commands.update_file(mock_cursor, mock_structure, MagicMock(), MagicMock())
        expected_calls = [
            call.execute("DELETE FROM main WHERE entry_id = '1A00'"),
            call.execute("INSERT INTO main VALUES(?, ?, ?)", ('1A00', 'data1', 'data2')),
            call.execute("DELETE FROM coils WHERE entry_id = '1A00'")
        ]
        
        mock_cursor.assert_has_calls(expected_calls)
    

