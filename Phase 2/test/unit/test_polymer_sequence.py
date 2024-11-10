"""
This script contains unit tests for testing methods in extract.py or polymer_sequence.py 
for helices, strands, and sheets.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, MagicMock
import gemmi 
from gemmi import cif

from polymer_sequence import Monomer, letter_code_3to1, sequence_3to1


def test_binary_search_target_found(mock_polymer_sequence):
    # sequence is defined in mock_polymer_sequence
    result = mock_polymer_sequence.binary_search(0, 4, 5)
    assert result == 2


def test_binary_search_target_found_at_start(mock_polymer_sequence):
    result = mock_polymer_sequence.binary_search(0, 3, 1)
    assert result == 0


def test_binary_search_target_found_at_end(mock_polymer_sequence):
    result = mock_polymer_sequence.binary_search(1, 4, 9)
    assert result == 4


def test_binary_search_target_not_found(mock_polymer_sequence):
    """ 
    Test that an Exception is raised when the target label cannot be found.
    """
    # target_label not in defined range
    with pytest.raises(Exception, match="Couldn't find index"):
        mock_polymer_sequence.binary_search(1, 4, 1)
    # target_label not in sequence
    with pytest.raises(Exception, match="Couldn't find index"):
        mock_polymer_sequence.binary_search(0, 4, 10)
        

def test_binary_search_invalid_range(mock_polymer_sequence):
    """ 
    Test that an Exception is raised when right index is less than left index.
    """
    with pytest.raises(Exception, match="Couldn't find index"):
        mock_polymer_sequence.binary_search(4, 0, 5)


def test_binary_search_empty_sequence(mock_polymer_sequence):
    """
    Test that an IndexError is raised when the sequence is empty.
    """
    mock_polymer_sequence.sequence = []

    with pytest.raises(IndexError):
        mock_polymer_sequence.binary_search(0, 4, 5)


def test_get_chain_subsequence_chain_not_in_start_indices(mock_polymer_sequence):
    """
    Test that an empty string is returned when the given chain is not found in chain_start_indices.
    """
    mock_chain_name = 'C'

    # start_id and end_id not required and can be mocked 
    result_sequence = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected_sequence = ''

    assert result_sequence == expected_sequence


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_larger_than_start(mock_binary_search, mock_polymer_sequence):
    """
    Test the get_chain_subsequence function when the end_index is 
    larger than or equal to the start_index.
    """
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [0, 10]
    
    # start_id and end_id not required and can be mocked 
    result = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected = ('ARNDCQEGHIX', 11)
 
    assert result == expected


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_equals_zero(mock_binary_search, mock_polymer_sequence):
    """
    Test the get_chain_subsequence function when the end_index is zero.
    """
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
    """
    Test the get_chain_subsequence function when the end_index is 
    less than the start_index and not zero.
    """ 
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [3, 1]

    # start_id and end_id not required and can be mocked 
    result = mock_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    
    # expected sequence is sliced backwards from start index
    expected = ('DNR', -3)

    assert result == expected


def test_get_chain_start_position(mock_polymer_sequence):
    mock_polymer_sequence.sequence = [Monomer('A', 1, 10, 'ALA', 'n')]
    result = mock_polymer_sequence.get_chain_start_position('A')
    assert result == 10


def test_get_chain_start_position_invalid_chain(mock_polymer_sequence):
    """ 
    Test that a KeyError is raised when the chain is not found in chain_start_indices.
    """
    mock_polymer_sequence.sequence = [Monomer('A', 1, 10, 'ALA', 'n')]
    with pytest.raises(KeyError, match='C'):
        mock_polymer_sequence.get_chain_start_position('C')


def test_get_chain_start_position_empty_sequence(mock_polymer_sequence):
    """ 
    Test that an IndexError is raised when the sequence is empty.
    """
    mock_polymer_sequence.sequence = []
    with pytest.raises(IndexError):
        mock_polymer_sequence.get_chain_start_position('A')


def test_get_chain_end_position(mock_polymer_sequence):
    mock_polymer_sequence.sequence = [Monomer('A', 1, 10, 'ALA', 'n'), Monomer('A', 1, 11, 'ARG', 'n')]
    mock_polymer_sequence.chain_end_indices = {'A': 1}
    result = mock_polymer_sequence.get_chain_end_position('A')
    assert result == 11


def test_get_chain_end_position_invalid_chain(mock_polymer_sequence):
    """ 
    Test that a KeyError is raised when the chain is not found in chain_end_indices
    """
    mock_polymer_sequence.sequence = [Monomer('A', 1, 10, 'ALA', 'n'), Monomer('A', 1, 11, 'ARG', 'n')]
    with pytest.raises(KeyError, match='C'):
        mock_polymer_sequence.get_chain_end_position('C')


def test_get_chain_end_position_empty_sequence(mock_polymer_sequence):
    """
    Test that an IndexError is raised when the sequence is empty.
    """
    mock_polymer_sequence.sequence = []
    with pytest.raises(IndexError):
        mock_polymer_sequence.get_chain_end_position('A')


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


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence_multiple_chains_error(mock_binary_search, mock_structure, mock_polymer_sequence, mock_helix):
    """
    Test the get_helix_sequence function when the start chain is not equal to the end_chain.
    """
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    helix_sequence = mock_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == "MULTIPLE CHAINS ERROR"


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence(mock_binary_search, mock_structure, mock_polymer_sequence, mock_strand):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = mock_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    expected = ('ARNDCQEGHIX', 11)
    
    assert result == expected


#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence_multiple_chains_error(mock_binary_search, mock_structure, mock_polymer_sequence, mock_strand):
    """
    Test the get_strand_sequence function when the start chain is not equal to the end_chain.
    """
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]

    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = mock_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    
    assert result == ("MULTIPLE CHAINS ERROR", -1)


def test_get_chain_sequence(mock_polymer_sequence):
    mock_chain_name = 'A'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = 'ARNDCQEGHIX'
    
    assert result_sequence == expected_sequence


def test_get_chain_sequence_chain_not_in_start_indices(mock_polymer_sequence):
    """
    Test that an empty string is returned when the given chain is not found in chain_start_indices.
    """
    mock_chain_name = 'C'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence


def test_get_chain_sequence_chain_not_in_end_indices(mock_polymer_sequence):
    """
    Test that an empty string is returned when the given chain is not found in chain_end_indices.
    """
    mock_chain_name = 'B'
    result_sequence = mock_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence


@patch.dict("polymer_sequence.three_to_one", {'AAA': 'A'})
def test_letter_code_3to1():
    known_polymer = 'AAA'    

    assert letter_code_3to1(known_polymer) == 'A'


@patch.dict("polymer_sequence.three_to_one", {'AAA': 'A'})
def test_letter_code_3to1_unknown_polymer():
    """
    Test that 'X' is returned when an unknown polymer is given. 
    """
    unknown_polymer = 'XXX'

    assert letter_code_3to1(unknown_polymer) == 'X'


def test_sequence_3to1():
    """
    Test the sequence_3to1 function for amino acid and DNA sequences.
    """
    amino_acid_seq = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'XXX']
    dna_seq = ['DA', 'DT', 'DG', 'DC']

    assert sequence_3to1(amino_acid_seq) == 'ARNDCQEGHIX'
    assert sequence_3to1(dna_seq) == 'ATGC'


def test_sequence_3to1_empty_sequence():
    """
    Test an empty string is returned when the given sequence is empty.
    """
    empty_sequence = []

    assert sequence_3to1(empty_sequence) == ''