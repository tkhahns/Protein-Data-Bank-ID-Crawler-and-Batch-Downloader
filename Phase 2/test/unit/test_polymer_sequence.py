"""
This script contains unit tests for testing methods in polymer_sequence.py.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""

import pytest
from unittest.mock import patch, MagicMock
import gemmi 
from gemmi import cif

from polymer_sequence import PolymerSequence, Monomer, letter_code_3to1, sequence_3to1, binary_search


def test_polymer_sequence_initialisation(mock_doc, fake_sequence_3to1):
    with patch("polymer_sequence.sequence_3to1", wraps=fake_sequence_3to1) as mock_sequence_3to1:
        # Initialize PolymerSequence with the mock CIF document
        polymer_sequence = PolymerSequence(mock_doc)

        # Last duplicate monomer is removed in the sequence
        expected_sequence = [
            Monomer("A", 1, 1, "ALA", "ALA", "n"),
            Monomer("A", 1, 2, "ARG", "?", "n"),
            Monomer("B", 1, 3, "ASN", "ASN", "y")
        ]
        
        assert polymer_sequence.sequence == expected_sequence
        assert polymer_sequence.bad_indices == [1]
        assert polymer_sequence.chain_start_indices == {"A": 0, "B": 2}
        assert polymer_sequence.chain_end_indices == {"A": 1, "B": 2}
        assert polymer_sequence.one_letter_code == "ARN"

        # Verify that duplicates for hetero entries were removed
        hetero_entries = [monomer for monomer in polymer_sequence.sequence if monomer.hetero == "y"]
        assert len(hetero_entries) == 1


def test_binary_search_target_found(test_polymer_sequence):
    # sequence defined in test_polymer_sequence
    result = test_polymer_sequence.binary_search(0, 4, 5)
    assert result == 4


def test_binary_search_target_found_at_start(test_polymer_sequence):
    result = test_polymer_sequence.binary_search(0, 3, 1)
    assert result == 0


def test_binary_search_target_found_at_end(test_polymer_sequence):
    result = test_polymer_sequence.binary_search(1, 10, 11)
    assert result == 10


def test_binary_search_target_not_found(test_polymer_sequence):
    """ 
    Test that an Exception is raised when the target label cannot be found.
    """
    # target_label not in defined range
    with pytest.raises(Exception, match="Couldn't find index"):
        test_polymer_sequence.binary_search(1, 4, 1)
    # target_label not in sequence
    with pytest.raises(Exception, match="Couldn't find index"):
        test_polymer_sequence.binary_search(0, 4, 10)
        

def test_binary_search_invalid_range(test_polymer_sequence):
    """ 
    Test that an Exception is raised when right index is less than left index.
    """
    with pytest.raises(Exception, match="Couldn't find index"):
        test_polymer_sequence.binary_search(4, 0, 5)


def test_binary_search_empty_sequence(test_polymer_sequence):
    """
    Test that an IndexError is raised when the sequence is empty.
    """
    test_polymer_sequence.sequence = []

    with pytest.raises(IndexError):
        test_polymer_sequence.binary_search(0, 4, 5)


def test_get_chain_subsequence_chain_not_in_start_indices(test_polymer_sequence):
    """
    Test that an empty string is returned when the given chain is not found in chain_start_indices.
    """
    mock_chain_name = 'C'

    # start_id and end_id not required and can be mocked 
    result_sequence = test_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected_sequence = ''

    assert result_sequence == expected_sequence


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_larger_than_start(mock_binary_search, test_polymer_sequence):
    """
    Test the get_chain_subsequence function when the end_index is 
    larger than or equal to the start_index.
    """
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [0, 10]
    
    # start_id and end_id not required and can be mocked 
    result = test_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    expected = ('ARNDCQEGHIX', 11)
 
    assert result == expected


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_equals_zero(mock_binary_search, test_polymer_sequence):
    """
    Test the get_chain_subsequence function when the end_index is zero.
    """
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [3, 0]

    # start_id and end_id not required and can be mocked 
    result = test_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    
    # expected sequence is sliced backwards from start index
    expected = ('DNRA', -4)

    assert result == expected


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_chain_subsequence_end_less_than_start_and_not_zero(mock_binary_search, test_polymer_sequence):
    """
    Test the get_chain_subsequence function when the end_index is 
    less than the start_index and not zero.
    """ 
    mock_chain_name = 'A'

    # mock return value of binary_search (start and end index)
    mock_binary_search.side_effect = [3, 1]

    # start_id and end_id not required and can be mocked 
    result = test_polymer_sequence.get_chain_subsequence(mock_chain_name, MagicMock(), MagicMock())
    
    # expected sequence is sliced backwards from start index
    expected = ('DNR', -3)

    assert result == expected


@patch("polymer_sequence.binary_search")
def test_get_chain_annotated_subsequence_end_larger_than_start_index(mock_binary_search, test_polymer_sequence, mock_span):
    """
    Test that the correct one letter sequence is returned when 
    the obtained end index is larger than or equal to the start index.
    """
    test_chain_string = 'ARNDCQEGH-IX'  # span_index_to_string_index will be [0, 1,..., 8, 10, 11]
    mock_binary_search.side_effect = [0, 9]  # mock start and end index 
    result = test_polymer_sequence.get_chain_annotated_subsequence(mock_span, test_chain_string, 1, 10)
    assert result == "ARNDCQEGH-I"


@patch("polymer_sequence.binary_search")
def test_get_chain_annotated_subsequence_empty_span(mock_binary_search, test_polymer_sequence):
    """
    Test that an empty string is returned when an empty span is given.
    """
    test_chain_string = 'ARNDCQEGH-IX'  # span_index_to_string_index will be [0, 1,..., 8, 10, 11]
    mock_binary_search.side_effect = [0, 9]  # mock start and end index 
    span = []
    result = test_polymer_sequence.get_chain_annotated_subsequence(span, test_chain_string, 1, 10)
    assert result == ""


@patch("polymer_sequence.binary_search")
def test_get_chain_annotated_subsequence_end_index_is_zero(mock_binary_search, test_polymer_sequence, mock_span):
    """
    Test that a reversed one letter sequence is returned when 
    the obtained end index equals zero.
    """
    test_chain_string = 'ARNDCQEGH-IX'  # span_index_to_string_index will be [0, 1,..., 8, 10, 11]
    mock_binary_search.side_effect = [3, 0]  # mock start and end index 
    result = test_polymer_sequence.get_chain_annotated_subsequence(mock_span, test_chain_string, 4, 1)
    assert result == "DNRA"  # sequence is reversed 


@patch("polymer_sequence.binary_search")
def test_get_chain_annotated_subsequence_end_smaller_than_start_index(mock_binary_search, test_polymer_sequence, mock_span):
    """
    Test that a reversed one letter sequence is returned when 
    the obtained end index is smaller than the start index but not equals zero.
    """
    test_chain_string = 'ARNDCQEGH-IX'  # span_index_to_string_index will be [0, 1,..., 8, 10, 11]
    mock_binary_search.side_effect = [3, 1]  # mock start and end index 
    result = test_polymer_sequence.get_chain_annotated_subsequence(mock_span, test_chain_string, 4, 2)
    assert result == "DNR"  # sequence is reversed 


def test_contains_unconfirmed_residues_residues_within_range(test_polymer_sequence):
    test_polymer_sequence.bad_indices = [5, 10]
    # start_id: 0, end_id: 11
    result = test_polymer_sequence.contains_unconfirmed_residues('A', 1, 11)
    assert result == 1 


def test_contains_unconfirmed_residues_no_residues_in_range(test_polymer_sequence):
    test_polymer_sequence.bad_indices = [5, 10]
    # start_id: 0, end_id: 5 
    result = test_polymer_sequence.contains_unconfirmed_residues('A', 1, 5)
    assert result == 0


def test_contains_unconfirmed_residues_no_bad_indices(test_polymer_sequence):
    test_polymer_sequence.bad_indices = []
    # start_id: 0, end_id: 11
    result = test_polymer_sequence.contains_unconfirmed_residues('A', 1, 11)
    assert result == 0


def test_get_chain_start_id(test_polymer_sequence):
    result = test_polymer_sequence.get_chain_start_id('A')
    assert result == 1


def test_get_chain_start_id_invalid_chain(test_polymer_sequence):
    """ 
    Test that a KeyError is raised when the chain is 
    not found in chain_start_indices.
    """
    with pytest.raises(KeyError, match='C'):
        test_polymer_sequence.get_chain_start_id('C')


def test_get_chain_start_id_empty_sequence(test_polymer_sequence):
    """ 
    Test that an IndexError is raised when the sequence is empty.
    """
    test_polymer_sequence.sequence = []
    with pytest.raises(IndexError):
        test_polymer_sequence.get_chain_start_id('A')


def test_get_chain_end_id(test_polymer_sequence):
    result = test_polymer_sequence.get_chain_end_id('A')
    assert result == 11


def test_get_chain_end_id_invalid_chain(test_polymer_sequence):
    """ 
    Test that a KeyError is raised when the chain is not found in chain_end_indices
    """
    with pytest.raises(KeyError, match='C'):
        test_polymer_sequence.get_chain_end_id('C')


def test_get_chain_end_id_empty_sequence(test_polymer_sequence):
    """
    Test that an IndexError is raised when the sequence is empty.
    """
    test_polymer_sequence.sequence = []
    with pytest.raises(IndexError):
        test_polymer_sequence.get_chain_end_id('A')


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence(mock_binary_search, mock_structure, test_polymer_sequence, mock_helix):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))

    expected_sequence = 'ARNDCQEGHIX'
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    helix_sequence = test_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == expected_sequence


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence_multiple_chains_error(mock_binary_search, mock_structure, test_polymer_sequence, mock_helix):
    """
    Test the get_helix_sequence function when the 
    start chain is not equal to the end_chain.
    """
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    helix_sequence = test_polymer_sequence.get_helix_sequence(mock_helix, mock_structure)
    
    assert helix_sequence == "MULTIPLE CHAINS ERROR"


@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence(mock_binary_search, mock_structure, test_polymer_sequence, mock_strand):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = test_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    expected = ('ARNDCQEGHIX', 11)
    
    assert result == expected


#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_strand_sequence_multiple_chains_error(mock_binary_search, mock_structure, test_polymer_sequence, mock_strand):
    """
    Test the get_strand_sequence function when the 
    start chain is not equal to the end_chain.
    """
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]

    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 10]

    result = test_polymer_sequence.get_strand_sequence(mock_strand, mock_structure)
    
    assert result == ("MULTIPLE CHAINS ERROR", 0)


def test_get_chain_sequence(test_polymer_sequence):
    mock_chain_name = 'A'
    result_sequence = test_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = 'ARNDCQEGHIX'
    
    assert result_sequence == expected_sequence


def test_get_chain_sequence_chain_not_in_start_indices(test_polymer_sequence):
    """
    Test that an empty string is returned when the given chain 
    is not found in chain_start_indices.
    """
    mock_chain_name = 'C'
    result_sequence = test_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence


def test_get_chain_sequence_chain_not_in_end_indices(test_polymer_sequence):
    """
    Test that an empty string is returned when the given chain 
    is not found in chain_end_indices.
    """
    mock_chain_name = 'B'
    result_sequence = test_polymer_sequence.get_chain_sequence(mock_chain_name)
    expected_sequence = ''

    assert result_sequence == expected_sequence


def test_span_binary_search_target_found(mock_span):
    """
    Test that the correct index is found when the target label is
    exactly in the middle of the given range. 
    """
    target_label = 5
    result = binary_search(mock_span, 0, len(mock_span) - 1, target_label)
    assert result == 4

def test_span_binary_search_recursion(mock_span):
    """
    Test that the correct index is found when a non-trivial 
    target label is given.
    """
    target_label = 8
    result = binary_search(mock_span, 0, len(mock_span) - 1, target_label)
    assert result == 7


def test_span_binary_search_target_below_range(mock_span):
    """
    Test that the given left index is returned when the target label
    is smaller than the smallest label_seq. 
    """
    target_label = 0
    result = binary_search(mock_span, 1, len(mock_span) - 1, target_label)
    assert result == 1


def test_span_binary_search_target_above_range(mock_span):
    """
    Test that the given right index is returned when the target label
    is larger than the largest label_seq. 
    """
    target_label = 11
    result = binary_search(mock_span, 0, len(mock_span) - 1, target_label)
    assert result == len(mock_span) - 1


def test_binary_search_invalid_range(mock_span):
    """
    Test that an Exception is raised when the given right index
    is smaller than the left index. 
    """
    with pytest.raises(Exception, match="Couldn't find index"):
        binary_search(mock_span, 7, 3, 5)
        

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