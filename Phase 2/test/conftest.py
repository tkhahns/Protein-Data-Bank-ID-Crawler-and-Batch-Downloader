import pytest
from unittest.mock import MagicMock, PropertyMock
import gemmi 
from gemmi import cif

from polymer_sequence import PolymerSequence

@pytest.fixture
def mock_structure():
    mock_structure = MagicMock(spec=gemmi.Structure)
    mock_structure.info = {'_entry.id': '1A00'}
    return mock_structure

@pytest.fixture
def mock_subchain():
    mock_subchain = MagicMock(spec=gemmi.ResidueSpan)
    mock_subchain.subchain_id.return_value = 'A'

    return mock_subchain

@pytest.fixture
def mock_polymer():
    mock_polymer = MagicMock()

    first_residue = MagicMock()
    last_residue = MagicMock()
    
    type(first_residue).label_seq = PropertyMock(return_value=1)
    type(last_residue).label_seq = PropertyMock(return_value=11)

    # mock_polymer.__getitem__.side_effect = [first_residue, last_residue] 
    mock_polymer.make_one_letter_sequence.return_value = 'ARNDCQEGHIX'
    mock_polymer.length.return_value = 11

    return mock_polymer

@pytest.fixture
def mock_chain(mock_subchain, mock_polymer):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'

    mock_chain.get_polymer.return_value = mock_polymer
    mock_chain.subchains.return_value = [mock_subchain, mock_subchain]

    return mock_chain

@pytest.fixture 
def mock_empty_polymer():
    mock_polymer = MagicMock(spec=gemmi.ResidueSpan)
    mock_polymer.make_one_letter_sequence.return_value = ''  
    mock_polymer.length.return_value = 0  

    return mock_polymer

@pytest.fixture
def mock_empty_chain(mock_subchain, mock_empty_polymer):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'

    mock_chain.get_polymer.return_value = mock_empty_polymer
    mock_chain.subchains.return_value = [mock_subchain]

    return mock_chain

@pytest.fixture
def mock_polymer_sequence():
    mock_doc = MagicMock(spec=cif.Document)
    polymer_sequence = PolymerSequence(mock_doc)

    polymer_sequence.one_letter_code = MagicMock()
    polymer_sequence.one_letter_code.__getitem__.return_value = 'ARNDCQEGHIX'

    polymer_sequence.chain_start_indices = {'A': 0, 'B': 0}
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
