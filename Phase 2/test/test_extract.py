"""
This script contains unit tests for testing methods in extract.py or polymer_sequence.py 
for complex types, helices, strands and sheets. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, MagicMock
import gemmi 
from gemmi import cif

import extract 
from extract import ComplexType
from polymer_sequence import PolymerSequence, sequence_3to1, letter_code_3to1
import mocks

@pytest.mark.skip(reason=None)
def test_insert_into_chain_table(mock_structure, mock_chain):
    """
        Tests the insert_into_chain_table function. 
    """
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_chain]

    mock_doc = MagicMock(spec=cif.Document)
    mock_polymer_sequence = MagicMock(spec=PolymerSequence)

    # mock unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = 'ARNDCQEGHIX'

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)

    expected = [
        ('1A00', 'A', 'A A', 'ARNDCQEGHIX', 'ARNDCQEGHIX', 1, 11, 11)
    ]

    assert result == expected

#@pytest.mark.skip(reason=None)
def test_insert_into_chain_table_start_end_pos_are_none(mock_structure, mock_empty_chain):
    """
        Tests that the insert_into_chain_table function returns None for start and end positions
        when the chain does not contain any polymers.
    """
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_empty_chain]
    
    mock_doc = MagicMock(spec=cif.Document)
    mock_polymer_sequence = MagicMock(spec=PolymerSequence)

    # mock empty unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = ''

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)

    expected = [
        ('1A00', 'A', 'A', '', '', None, None, 0)
    ]

    assert result == expected

#@pytest.mark.skip(reason=None)
def test_get_chain_sequence(mock_polymer_sequence):
    mock_chain_name = 'A'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = 'ARNDCQEGHIX'
    
    # check that one_letter_code was indexed once
    mock_polymer_sequence.one_letter_code.__getitem__.assert_called_once()
    assert result_sequence == expected_sequence

#@pytest.mark.skip(reason=None)
def test_get_chain_sequence_chain_not_in_start_indices(mock_polymer_sequence):
    mock_chain_name = 'C'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence

@pytest.mark.xfail(reason="condition checking for chain in end_indices absent") 
def test_get_chain_sequence_chain_not_in_end_indices(mock_polymer_sequence):
    mock_chain_name = 'B'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence

@patch.dict("polymer_sequence.three_to_one", {'AAA': 'A'})
def test_letter_code_3to1():
    known_polymer = 'AAA'
    unknown_polymer = 'XXX'
    
    assert letter_code_3to1(known_polymer) == 'A'
    assert letter_code_3to1(unknown_polymer) == 'X'

@pytest.mark.xfail(reason="last monomer is omitted") 
def test_sequence_3to1():
    amino_acid_sequence = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'XXX']
    result = sequence_3to1(amino_acid_sequence)
    expected = 'ARNDCQEGHIX'

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
    mock_binary_search.side_effect = [0, 9]

    helix_sequence = mock_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == "MULTIPLE CHAINS ERROR"

#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence(mock_binary_search, mock_structure, mock_polymer_sequence, mock_strand):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))

    expected = ('ARNDCQEGHIX', 10)
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 9]

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
    mock_binary_search.side_effect = [0, 9]

    result = mock_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    
    assert result == ("MULTIPLE CHAINS ERROR", -1)

#@pytest.mark.skip(reason=None)
@pytest.mark.parametrize("complex_type", ComplexType)
def test_get_complex_type(complex_type):
    result = extract.get_complex_type(mocks.MockComplexType(complex_type).get_struct())
    expected = complex_type 
    assert(result == expected)

#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence', autospec = True)
def test_sheet_table(mock_polymer_sequence):
    mock_sheet_table = mocks.MockSheetTable()
    mock_struct, mock_doc = mock_sheet_table.get_doc()

    expected = [
        ("0A0A", "A", 3, "PPA"), 
        ("0A0A", "B", 3, "PA")
    ]
    result = extract.insert_into_sheet_table(mock_struct, mock_doc, mock_polymer_sequence)

    assert (result == expected)

