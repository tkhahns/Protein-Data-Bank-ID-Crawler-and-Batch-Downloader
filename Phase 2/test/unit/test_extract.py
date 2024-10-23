"""
This script contains unit tests for testing methods in extract.py or polymer_sequence.py 
for complex types and chains. 
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
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


@pytest.mark.parametrize("complex_type", extract.ComplexType)
def test_get_complex_type(complex_type, mock_structure, mock_entities):
    """
    Test that the function returns the correct complex type based on 
    the types of the entities.
    """
    mock_structure.entities = mock_entities(complex_type=complex_type)

    result = extract.get_complex_type(mock_structure)
    expected = complex_type 

    assert(result == expected)


@patch('extract.get_complex_type')
def test_insert_into_main_table(mock_complex_type, mock_structure, mock_doc, mock_chain):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_block = MagicMock(spec=cif.Block)
    mock_block.find_value.return_value = "'mock_org'"
    mock_doc.sole_block.return_value = mock_block
    mock_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'NucleicAcid', 'mock_title', 'mock_org', 'A A', 'P 1',
         1, 1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ]

    assert result == expected

@patch('extract.get_complex_type')
def test_insert_into_main_table_source_org_is_none(mock_get_complex_type, mock_structure, mock_doc, mock_chain):
    """
    Test that source organism is an empty string when its value is None.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_block = MagicMock(spec=cif.Block)
    mock_block.find_value.return_value = None
    mock_doc.sole_block.return_value = mock_block
    mock_get_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'NucleicAcid', 'mock_title', '', 'A A', 'P 1',
         1, 1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ]

    assert result == expected


@patch('extract.get_complex_type')
def test_insert_into_main_table_z_value_absent(mock_get_complex_type, mock_structure, mock_doc, mock_chain):
    """
    Test that z value is an empty string when the attribute is absent.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_block = MagicMock(spec=cif.Block)
    mock_doc.sole_block.return_value = mock_block
    mock_block.find_value.return_value = "'mock_org'"
    mock_get_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.info = {'_entry.id': '1A00', '_struct.title': 'mock_title'}  # z_value removed
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'NucleicAcid', 'mock_title', 'mock_org', 'A A', 'P 1',
         '', 1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ]

    assert result == expected


def test_insert_into_main_table_invalid_block(mock_structure):
    """
    Test that an error is raised when the cif Document contains 
    no blocks or more than one block.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    doc_one = cif.Document()  # has two blocks
    doc_one.add_new_block('block_one')
    doc_one.add_new_block('block_two')
    doc_two = cif.Document()  # one block 

    with pytest.raises(RuntimeError, match="single data block expected, got 2"):
        extract.insert_into_main_table(mock_structure, doc_one, mock_polymer_sequence)
    
    with pytest.raises(IndexError):
        extract.insert_into_main_table(mock_structure, doc_two, mock_polymer_sequence)


def test_insert_into_experimental_table(mock_structure, mock_doc):
    mock_values = {
        '_exptl_crystal.density_Matthews': 1.0, 
        '_exptl_crystal.density_percent_sol': 1.0,
        '_exptl_crystal_grow.method': 'mock_growth_method', 
        '_exptl_crystal_grow.pdbx_details': 'mock_growth_proc',
        '_exptl_crystal_grow.apparatus': 'mock_growth_apparatus', 
        '_exptl_crystal_grow.atmosphere': 'mock_growth_atmosophere', 
        '_exptl_crystal_grow.pH': 7.0, 
        '_exptl_crystal_grow.temp': 200.0
    }

    def side_effect(arg):
        return mock_values[arg]
    
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_block = MagicMock(spec=cif.Block)
    mock_doc.sole_block.return_value = mock_block
    mock_block.find_value.side_effect = side_effect

    result = extract.insert_into_experimental_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 1.0, 1.0, 'mock_growth_method', 'mock_growth_proc', 'mock_growth_apparatus', 
         'mock_growth_atmosophere', 7.0, 200.0)
    ]

    assert result == expected


def test_insert_into_experimental_table_missing_data(mock_structure, mock_doc):
    """
    Test that each data item is an empty string when its value is None.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_block = MagicMock(spec=cif.Block)
    mock_doc.sole_block.return_value = mock_block
    mock_block.find_value.return_value = None

    result = extract.insert_into_experimental_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', '', '', '', '', '', '', '', '')
    ]

    assert result == expected


def test_insert_into_experimental_table_invalid_block(mock_structure):
    """
    Test that an error is raised when the cif Document contains 
    no blocks or more than one block.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    doc_one = cif.Document()
    doc_one.add_new_block('block_one')
    doc_one.add_new_block('block_two')
    doc_two = cif.Document()

    with pytest.raises(RuntimeError, match="single data block expected, got 2") as errinfo_one:
        extract.insert_into_experimental_table(mock_structure, doc_one, mock_polymer_sequence)
    
    with pytest.raises(IndexError):
        extract.insert_into_experimental_table(mock_structure, doc_two, mock_polymer_sequence)
    

def test_insert_into_entity_table(mock_structure, mock_doc, mock_entity):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_block = MagicMock(spec=cif.Block)
    mock_doc.sole_block.return_value = mock_block
    mock_block.find_loop.return_value = ["mock_entity_one", "'mock_entity_two'", None]

    mock_entity_one = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideD, subchains=['A'])
    mock_entity_two = mock_entity(gemmi.EntityType.NonPolymer) # unknown polymer type, no subchains
    mock_entity_three = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideD, subchains=['A', 'B'])
    mock_structure.entities = [mock_entity_one, mock_entity_two, mock_entity_three]

    result = extract.insert_into_entity_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', '1', 'mock_entity_one', 'Polymer', 'PeptideD', 'A'), 
        ('1A00', '1', 'mock_entity_two', 'NonPolymer', 'Unknown', ''),
        ('1A00', '1', '', 'Polymer', 'PeptideD', 'A B'),
    ]

    assert result == expected
   

def test_insert_into_entity_table_loop_not_found(mock_structure, mock_doc, mock_entity):
    """
    Test that a list containing one tuple representing a single entity
    is returned when a loop is not found.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_block = MagicMock(spec=cif.Block)
    mock_doc.sole_block.return_value = mock_block
    mock_block.find_loop.return_value = []
    mock_block.find_value.return_value = "'mock_entity'"

    mock_entity = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideD, subchains=['A', 'B'])
    mock_structure.entities = [mock_entity]

    result = extract.insert_into_entity_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', '1', 'mock_entity', 'Polymer', 'PeptideD', 'A B')
    ]

    assert result == expected

 
def test_insert_into_entity_table_invalid_block(mock_structure):
    """
    Test that an error is raised when the cif Document contains 
    no blocks or more than one block.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    doc_one = cif.Document()
    doc_one.add_new_block('block_one')
    doc_one.add_new_block('block_two')
    doc_two = cif.Document()

    with pytest.raises(RuntimeError, match="single data block expected, got 2"):
        extract.insert_into_entity_table(mock_structure, doc_one, mock_polymer_sequence)
    
    with pytest.raises(IndexError):
        extract.insert_into_entity_table(mock_structure, doc_two, mock_polymer_sequence)


# def test_insert_into_subchain_table():
#     pass


def test_insert_into_chain_table(mock_structure, mock_doc, mock_chain):
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_chain]

    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_polymer_sequence.get_chain_start_position.return_value = 1
    mock_polymer_sequence.get_chain_end_position.return_value = 11

    # mock unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = 'ARNDCQEGHIX'

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'A', 'A A', 'ARNDCQEGHIX', 'ARNDCQEGHIX', 1, 11, 11)
    ]

    assert result == expected


def test_insert_into_chain_table_start_end_pos_are_none(mock_structure, mock_doc, mock_empty_chain):
    """
    Test that start and end positions are None when the chain does not contain any polymers.
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


def test_insert_into_helix_table_different_start_and_end_chain(mock_structure, mock_doc, mock_helix):    
    """
    Test that chain names are concatenated when the start and end chain are different. 
    """
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


# def test_insert_into_coil_table():
#     pass
