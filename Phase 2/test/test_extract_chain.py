"""
This script contains unit tests for testing methods in extract.py or polymer_sequence.py 
for complex types and chains. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, MagicMock
from gemmi import cif

import extract 
import polymer_sequence


@pytest.mark.parametrize("complex_type", extract.ComplexType)
def test_get_complex_type(complex_type, mock_structure, mock_entities):
    mock_structure.entities = mock_entities(complex_type=complex_type)
    result = extract.get_complex_type(mock_structure)
    expected = complex_type 

    assert(result == expected)


def test_insert_into_chain_table(mock_structure, mock_doc, mock_chain):
    """
        Tests the insert_into_chain_table function. 
    """
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_chain]

    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    # create mock dictionaries for start and end indices
    mock_polymer_sequence.chain_start_indices = MagicMock()
    mock_polymer_sequence.chain_end_indices = MagicMock()

    # set up return values for mock dictionaries 
    mock_polymer_sequence.chain_start_indices.__getitem__.side_effect = lambda key: {'A': 0}[key]
    mock_polymer_sequence.chain_end_indices.__getitem__.side_effect = lambda key: {'A': 10}[key]

    # mock unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = 'ARNDCQEGHIX'

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)

    expected = [
        ('1A00', 'A', 'A A', 'ARNDCQEGHIX', 'ARNDCQEGHIX', 0, 10, 11)
    ]

    # check if start and end indices dictionaries are accessed with correct keys
    mock_polymer_sequence.chain_start_indices.__getitem__.assert_called_once_with('A')
    mock_polymer_sequence.chain_end_indices.__getitem__.assert_called_once_with('A')

    assert result == expected


def test_insert_into_chain_table_start_end_pos_are_none(mock_structure, mock_doc, mock_empty_chain):
    """
        Tests that the insert_into_chain_table function returns None for start and end positions
        when the chain does not contain any polymers.
    """
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_empty_chain]
    
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    # mock empty unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = ''

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)

    expected = [
        ('1A00', 'A', 'A', '', '', None, None, 0)
    ]

    assert result == expected


def test_get_chain_sequence(mock_polymer_sequence):
    mock_chain_name = 'A'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = 'ARNDCQEGHIX'
    
    assert result_sequence == expected_sequence


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
    
    assert polymer_sequence.letter_code_3to1(known_polymer) == 'A'
    assert polymer_sequence.letter_code_3to1(unknown_polymer) == 'X'


@pytest.mark.xfail(reason="last monomer is omitted") 
def test_sequence_3to1():
    amino_acid_seq = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'XXX']
    result_aa_seq = polymer_sequence.sequence_3to1(amino_acid_seq)
    expected_aa_seq = 'ARNDCQEGHIX'

    dna_seq = ['DA', 'DT', 'DG', 'DC']
    result_dna_seq = polymer_sequence.sequence_3to1(dna_seq)
    expected_dna_seq = 'ATGC'
    
    assert result_aa_seq == expected_aa_seq
    assert result_dna_seq == expected_dna_seq

