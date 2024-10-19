"""
This script contains unit tests for testing methods in extract.py or polymer_sequence.py 
for helices, strands, and sheets.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, MagicMock, call
import gemmi 
from gemmi import cif

import extract
import polymer_sequence


def test_sense_sequence(mock_sheet):
    mock_strand_1 = MagicMock(spec=gemmi.Sheet.Strand, sense = 1)
    mock_strand_2 = MagicMock(spec=gemmi.Sheet.Strand, sense = -1)
    mock_strand_3 = MagicMock(spec=gemmi.Sheet.Strand, sense = 0)

    mock_sheet.strands = [mock_strand_1, mock_strand_2, mock_strand_3]

    result = extract.sense_sequence(mock_sheet)

    expected = 'PA'

    assert result == expected 


def test_get_chain_subsequence_chain_not_in_start_indices(mock_polymer_sequence):
    mock_chain_name = 'C'

    # start_id and end_id not required and can be mocked 
    result_sequence = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected_sequence = ''

    assert result_sequence == expected_sequence


@patch('polymer_sequence.PolymerSequence.binary_search')
#@pytest.mark.skip(reason=None)
def test_get_chain_subsequence_end_larger_than_start(mock_binary_search, mock_polymer_sequence):
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [0, 10]
    
    # start_id and end_id not required and can be mocked 
    result = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected = ('ARNDCQEGHIX', 11)
 
    assert result == expected


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_equals_zero(mock_binary_search, mock_polymer_sequence):
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [3, 0]

    # start_id and end_id not required and can be mocked 
    result = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    
    # expected sequence is sliced backwards from start index
    expected = ('DNRA', -4)

    assert result == expected


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_less_than_start_and_not_zero(mock_binary_search, mock_polymer_sequence):
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [3, 1]

    # start_id and end_id not required and can be mocked 
    result = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    
    # expected sequence is sliced backwards from start index
    expected = ('DNR', -3)

    assert result == expected


def test_insert_into_helix_table_different_start_and_end_chain(mock_structure, mock_doc, mock_helix):    
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    # sequence_code
    mock_polymer_sequence.get_helix_sequence.return_value = 'ARNDCQEGHIX'
    mock_structure.helices = [mock_helix]

    # different chain and end_chain
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # start_pos, end_pos 
    mock_start_chain.__getitem__.return_value = [MagicMock(label_seq=1)]
    mock_end_chain.__getitem__.return_value = [MagicMock(label_seq=11)]

    result = extract.insert_into_helix_table(mock_structure, mock_doc, mock_polymer_sequence)
    
    expected = [
        ('1A00', 1, 'A B', 'ARNDCQEGHIX', 1, 11, 11)
    ]
    
    # check if find_cra was called with helix start and end 
    expected_calls = [call.find_cra(mock_helix.start), call.find_cra(mock_helix.end)]
    mock_structure[0].assert_has_calls(expected_calls)
 
    # check if mock chains were accessed with helix auth labels
    mock_start_chain.__getitem__.assert_called_once_with('1')
    mock_end_chain.__getitem__.assert_called_once_with('11')

    assert result == expected


#@pytest.mark.skip(reason=None)
def test_insert_into_helix_table_same_start_and_end_chain(mock_structure, mock_doc, mock_helix):    
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    # sequence_code
    mock_polymer_sequence.get_helix_sequence.return_value = 'ARNDCQEGHIX'

    mock_structure.helices = [mock_helix]

    # same chain and end_chain
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_chain), MagicMock(chain=mock_chain)]
    
    # start_pos, end_pos 
    mock_chain.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=11)]]

    result = extract.insert_into_helix_table(mock_structure, mock_doc, mock_polymer_sequence)
    
    expected = [
        ('1A00', 1, 'A', 'ARNDCQEGHIX', 1, 11, 11)
    ]
    
    # check if mock chain was accessed with helix auth labels
    expected_calls = [call.__getitem__('1'), call.__getitem__('11')]
    mock_chain.assert_has_calls(expected_calls)

    assert result == expected
    

#@pytest.mark.skip(reason=None)
@patch('extract.sense_sequence')
def test_insert_into_sheet_table(mock_sense_sequence, mock_structure, mock_doc, mock_sheet, mock_strand):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_structure.sheets = [mock_sheet]
    mock_sheet.strands = [mock_strand]

    mock_sense_sequence.return_value = 'P'

    result = extract.insert_into_sheet_table(mock_structure, mock_doc, mock_polymer_sequence)

    expected = [
        ('1A00', 'A', 1, 'P')
    ]

    assert result == expected


#@pytest.mark.skip(reason=None)
def test_insert_into_strand_table(mock_structure, mock_doc, mock_sheet, mock_strand):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_structure.sheets = [mock_sheet]
    mock_sheet.strands = [mock_strand]

    # sequence_code, length
    mock_polymer_sequence.get_strand_sequence.return_value = ['ARNDCQEGHIX', 11]
    
    # chain and end_chain
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # start_pos, end_pos 
    mock_start_chain.__getitem__.return_value = [MagicMock(label_seq=1)]
    mock_end_chain.__getitem__.return_value = [MagicMock(label_seq=11)]

    result = extract.insert_into_strand_table(mock_structure, mock_doc, mock_polymer_sequence)
    
    expected = [
        ('1A00', 'A', '1', 'A', 'ARNDCQEGHIX', 1, 11, 11)
    ]
    
    # check if find_cra was called with strand start and end 
    expected_calls = [call.find_cra(mock_strand.start), call.find_cra(mock_strand.end)]
    mock_structure[0].assert_has_calls(expected_calls)

    # check if mock chains were accessed with strand auth labels
    mock_start_chain.__getitem__.assert_called_once_with('1')
    mock_end_chain.__getitem__.assert_called_once_with('11')

    assert result == expected

 
#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence(mock_binary_search, mock_structure, mock_polymer_sequence, mock_helix):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))

    expected_sequence = 'ARNDCQEGHIX'
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    helix_sequence = mock_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == expected_sequence


#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence_multiple_chains_error(mock_binary_search, mock_structure, mock_polymer_sequence, mock_helix):
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    helix_sequence = mock_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == "MULTIPLE CHAINS ERROR"

#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence(mock_binary_search, mock_structure, mock_polymer_sequence, mock_strand):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))

    expected = ('ARNDCQEGHIX', 11)
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = mock_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    
    assert result == expected


#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence_multiple_chains_error(mock_binary_search, mock_structure, mock_polymer_sequence, mock_strand):
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]

    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = mock_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    
    assert result == ("MULTIPLE CHAINS ERROR", -1)

