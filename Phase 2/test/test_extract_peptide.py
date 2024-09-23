import pytest
from unittest.mock import patch, MagicMock
import gemmi 

import extract
import mocks

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
@patch('polymer_sequence.PolymerSequence', autospec = True)
def test_insert_into_sheet_table(mock_polymer_sequence):
    mock_sheet_table = mocks.MockSheetTable()
    mock_struct, mock_doc = mock_sheet_table.get_doc()

    expected = [
        ("0A0A", "A", 3, "PPA"), 
        ("0A0A", "B", 3, "PA")
    ]
    result = extract.insert_into_sheet_table(mock_struct, mock_doc, mock_polymer_sequence)

    assert (result == expected)
