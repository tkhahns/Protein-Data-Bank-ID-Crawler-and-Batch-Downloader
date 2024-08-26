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
from polymer_sequence import PolymerSequence
import mocks

@pytest.fixture
def mock_structure():
    mock_structure = MagicMock(spec=gemmi.Structure)

    return mock_structure

@pytest.fixture
def mock_polymer_sequence():
    mock_doc = MagicMock(spec=cif.Document)
    polymer_sequence = PolymerSequence(mock_doc)
    polymer_sequence.one_letter_code = 'ADQCMVIRQMPWPCATWRPRNEYLMHKHRRGMDPCYIGPFPGCQIGINIK'

    polymer_sequence.chain_start_indices = {'A': 0}
    polymer_sequence.chain_end_indices = {'A': len(polymer_sequence.one_letter_code) - 1}

    return polymer_sequence

@pytest.fixture
def mock_helix():
    mock_helix = MagicMock(spec=gemmi.Helix)
    mock_helix.start.res_id.seqid.num = 1
    mock_helix.start.res_id.seqid.icode = ''
    mock_helix.end.res_id.seqid.num = 10
    mock_helix.end.res_id.seqid.icode = ''

    return mock_helix

@pytest.fixture
def mock_strand():
    mock_strand = MagicMock(spec=gemmi.Sheet.Strand)
    mock_strand.start.res_id.seqid.num = 41
    mock_strand.start.res_id.seqid.icode = ''
    mock_strand.end.res_id.seqid.num = 50
    mock_strand.end.res_id.seqid.icode = ''

    return mock_strand

#@pytest.mark.skip(reason=None)
@patch('polymer_sequence.PolymerSequence.binary_search')
def test_get_helix_sequence(mock_binary_search, mock_structure, mock_polymer_sequence, mock_helix):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.return_value.chain = mock_chain
    
    mock_chain.__getitem__.side_effect = lambda x: MagicMock(label_seq=int(x))

    expected_sequence = "ADQCMVIRQM"
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [0, 9]

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

    expected = ("PGCQIGINIK", 10)
    
    # mock return values of binary_search function to chain_start and chain_end
    mock_binary_search.side_effect = [40, 49]

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

