"""
This script contains unit tests for testing methods in extract.py.
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
    mock_block.find_value.side_effect = ["'mock_org'", "2000-01-01"]
    mock_doc.sole_block.return_value = mock_block
    mock_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ("1A00", "NucleicAcid", "mock_title", "mock_org", "2000-01-01", "A A", "P 1",
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
    mock_block.find_value.side_effect = [None, "2000-01-01"]
    mock_doc.sole_block.return_value = mock_block
    mock_get_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ("1A00", "NucleicAcid", "mock_title", "", "2000-01-01", "A A", "P 1",
         1, 1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ]

    assert result == expected

@patch('extract.get_complex_type')
def test_insert_into_main_table_revision_date_is_none(mock_get_complex_type, mock_structure, mock_doc, mock_chain):
    """
    Test that source organism is an empty string when its value is None.
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_block = MagicMock(spec=cif.Block)
    mock_block.find_value.side_effect = ["'mock_org'", None]
    mock_block.find_loop.return_value = ["2000-01-01"]
    mock_doc.sole_block.return_value = mock_block
    mock_get_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ("1A00", "NucleicAcid", "mock_title", "mock_org", "2000-01-01", "A A", "P 1",
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
    mock_block.find_value.side_effect = ["'mock_org'", "2000-01-01"]
    mock_get_complex_type.return_value = extract.ComplexType.NucleicAcid
    mock_structure.info = {'_entry.id': '1A00', '_struct.title': 'mock_title'}  # z_value removed
    mock_structure.__getitem__.return_value = [mock_chain, mock_chain]
    
    result = extract.insert_into_main_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ("1A00", "NucleicAcid", "mock_title", "mock_org", "2000-01-01", "A A", "P 1",
         "", 1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
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


def test_insert_into_subchain_table(mock_structure, mock_doc, mock_entity, mock_subchain, mock_chain):
    def mock_get_item(index):
        if index == 0:
            return mock_first_residue
        elif index == -1:
            return mock_last_residue
        
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_entity = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideD, ['A'])
    mock_structure.entities = [mock_entity]

    # mock subchain length to a non-zero value
    mock_subchain.__len__.return_value = 11  
    mock_structure[0].get_subchain.return_value = mock_subchain
    mock_structure[0].get_parent_of.return_value = mock_chain

    # mock start and end positions
    mock_first_residue = MagicMock(label_seq=1)
    mock_last_residue = MagicMock(label_seq=11)
    mock_subchain.__getitem__.side_effect = mock_get_item 

    # mock unannotated sequence and unconfirmed 
    mock_polymer_sequence.get_chain_subsequence.return_value = ['ARNDCQEGHIX', 11]

    result = extract.insert_into_subchain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', '1', 'A', 'A', 'ARNDCQEGHIX', 'ARNDCQEGHIX', 1, 11, 11)
    ]
    
    assert result == expected
  

def test_insert_into_subchain_table_invalid_entity(mock_structure, mock_doc, mock_entity):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_invalid_entity = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.SaccharideD, ['A'])  # skipped in insert method
    mock_structure.entities = [mock_invalid_entity]

    result = extract.insert_into_subchain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = []

    assert result == expected


def test_insert_into_subchain_table_subchain_length_is_zero(mock_structure, mock_doc, mock_entity, mock_subchain, mock_chain):
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    mock_entity = mock_entity(gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideD, ['A'])
    mock_structure.entities = [mock_entity]

    # mock subchain length to zero 
    mock_subchain.__len__.return_value = 0
    mock_structure[0].get_subchain.return_value = mock_subchain
    mock_structure[0].get_parent_of.return_value = mock_chain

    result = extract.insert_into_subchain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = []

    assert result == expected


def test_insert_into_chain_table(mock_structure, mock_doc, mock_chain):
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_chain]

    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    mock_polymer_sequence.get_chain_start_id.return_value = 1
    mock_polymer_sequence.get_chain_end_id.return_value = 11

    # mock unannotated sequence and unconfirmed 
    mock_polymer_sequence.get_chain_sequence.return_value = 'ARNDCQEGHIX'
    mock_polymer_sequence.contains_unconfirmed_residues.return_value = 0

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'A', 'A A', 0, 'ARNDCQEGHIX', 'ARNDCQEGHIX', 1, 11, 11)
    ]

    assert result == expected


def test_insert_into_chain_table_empty_polymer(mock_structure, mock_doc, mock_empty_chain):
    """
    Test that start position, end position and unconfirmed are None 
    when the chain contains an empty polymer with no residues. 
    """
    # assign chain to mock_structure
    mock_structure.__getitem__.return_value = [mock_empty_chain]
    
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)
    
    # mock empty unannotated sequence
    mock_polymer_sequence.get_chain_sequence.return_value = ''

    result = extract.insert_into_chain_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 'A', 'A', None, '', '', None, None, 0)
    ]

    assert result == expected


def test_insert_into_helix_table_different_start_and_end_chain(mock_structure, mock_doc, mock_helix):    
    """
    Test that chain names are concatenated when the start and end chain are different. 
    """
    mock_polymer_sequence = MagicMock(spec=polymer_sequence.PolymerSequence)

    # mock helix sequence 
    mock_polymer_sequence.get_helix_sequence.return_value = 'ARNDCQEGHIX'
    mock_structure.helices = [mock_helix]

    # different chain and end_chain
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # mock start and end positions
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

    # mock helix sequence 
    mock_polymer_sequence.get_helix_sequence.return_value = 'ARNDCQEGHIX'
    mock_structure.helices = [mock_helix]

    # same chain and end_chain
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_chain), MagicMock(chain=mock_chain)]
    
    # mock start and end positions
    mock_chain.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=11)]]

    result = extract.insert_into_helix_table(mock_structure, mock_doc, mock_polymer_sequence)
    expected = [
        ('1A00', 1, 'A', 'ARNDCQEGHIX', 1, 11, 11)
    ]
    
    # check if find_cra was called with helix start and end 
    expected_calls = [call.find_cra(mock_helix.start), call.find_cra(mock_helix.end)]
    mock_structure[0].assert_has_calls(expected_calls)

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

    # mock strand sequence and length 
    mock_polymer_sequence.get_strand_sequence.return_value = ['ARNDCQEGHIX', 11]
    
    # chain and end_chain
    mock_start_chain = MagicMock(spec=gemmi.Chain)
    mock_start_chain.name = 'A'
    mock_end_chain = MagicMock(spec=gemmi.Chain)
    mock_end_chain.name = 'B'
    mock_structure[0].find_cra.side_effect = [MagicMock(chain=mock_start_chain), MagicMock(chain=mock_end_chain)]
    
    # mock start and end positions
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

class TestCoilExtractor:
    @pytest.fixture
    @staticmethod
    def mock_helix():
        helix = MagicMock()
        helix.start.res_id.seqid.num = 1
        helix.start.res_id.seqid.icode = ''
        helix.end.res_id.seqid.num = 5
        helix.end.res_id.seqid.icode = ''
        return helix
    
    @pytest.fixture
    @staticmethod
    def mock_sheet():
        strand = MagicMock()
        strand.name = "1"
        strand.start.res_id.seqid.num = 6
        strand.start.res_id.seqid.icode = ''
        strand.end.res_id.seqid.num = 7
        strand.end.res_id.seqid.icode = ''
        sheet = MagicMock()
        sheet.name = "A"
        sheet.strands = [strand]
        return sheet
    
    @pytest.fixture
    @staticmethod
    def mock_chain_a():
        chain = MagicMock()
        chain.name = "A"
        chain.get_polymer.return_value.make_one_letter_sequence.return_value = "ARNDCQEGHIX"
        chain.get_polymer.return_value.first_conformer.return_value = ["res1", "res2"]
        return chain
    
    @pytest.fixture
    @staticmethod
    def mock_chain_b():
        chain = MagicMock()
        chain.name = "B"
        return chain

    @pytest.fixture
    @staticmethod
    def mock_structure(mock_helix, mock_sheet, mock_chain):
        structure = MagicMock(spec=gemmi.Structure)
        structure.info = {"_entry.id": "1A00"}
        structure.helices = [mock_helix]
        structure.sheets = [mock_sheet]
        structure[0].find_chain.side_effect = lambda x: mock_chain if x == "A" else None
        return structure

    @pytest.fixture 
    @staticmethod 
    def mock_polymer_sequence():
        polymer_sequence = MagicMock()
        polymer_sequence.one_letter_code = 'ARNDCQEGHIX'

        polymer_sequence.chain_start_indices = {"A": 1, "B": 1}
        polymer_sequence.chain_end_indices = {"A": 10, "B": 8}

        polymer_sequence.get_chain_start_id.side_effect = lambda chain: 1 if chain == "A" else 1
        polymer_sequence.get_chain_end_id.side_effect = lambda chain: 10 if chain == "A" else 8

        polymer_sequence.get_chain_subsequence.return_value = ["SUBSEQ", 6]
        polymer_sequence.get_chain_annotated_subsequence.return_value = "SUBSEQ"
        polymer_sequence.contains_unconfirmed_residues.return_value = 0

        return polymer_sequence

    def test_insert_into_coil_table(self, mock_structure, mock_polymer_sequence, mock_chain_a, mock_chain_b):
        """
        Test valid extraction of coils with helices and sheets.
        """
        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=5)]]
        mock_chain_b.__getitem__.side_effect = [[MagicMock(label_seq=6)], [MagicMock(label_seq=7)]]  

        # mock start and end chains -- here, the helix and sheet are each in one chain
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
        ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "A", 0, "SUBSEQ", "SUBSEQ", 6, 10, 6), 
            ("1A00", 2, "B", 0, "SUBSEQ", "SUBSEQ", 1, 5, 6),
            ("1A00", 3, "B", 0, "SUBSEQ", "SUBSEQ", 8, 8, 6)
        ]

        # Expected output: 3 coils (before and after helix/sheet)
        assert len(result) == 3
        assert result == expected

    def test_insert_into_coil_table_ill_defined_helix(self, mock_structure, mock_polymer_sequence, mock_chain_a, mock_chain_b, capsys):
        """
        Test that no coils are extracted when a helix spans multiple chains
        (i.e. is ill-defined).
        """
        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)]]
        mock_chain_b.__getitem__.return_value = [MagicMock(label_seq=5), [MagicMock(label_seq=6)], [MagicMock(label_seq=7)]]

        # mock start and end chains -- here, the helix spans across two chains 
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
            ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = []
        expected_output =  "Helix 0 in protein 1A00 is ill-defined. Unable to extract random coils."

        # Expected output: 0 coils extracted
        assert result == expected
        # Expected print output 
        assert expected_output in capsys.readouterr().out 

    
    def test_insert_into_coil_table_ill_defined_strand(self, mock_structure, mock_polymer_sequence, mock_chain_a, mock_chain_b, capsys):
        """
        Test that no coils are extracted when a strand spans multiple chains
        (i.e. is ill-defined).
        """
        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=5)],
                                              [MagicMock(label_seq=6)]]
        mock_chain_b.__getitem__.return_value = [MagicMock(label_seq=7)]

        # mock start and end chains -- here, the sheet spans across two chains 
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b)
            ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = []
        expected_output =  "Strand 1 in sheet A in protein 1A00 is ill-defined. Unable to extract random coils."

        # Expected output: 0 coils extracted
        assert result == expected
        # Expected print output 
        assert expected_output in capsys.readouterr().out 

    
    def test_insert_into_coil_table_no_secondary_structures(self, mock_structure, mock_polymer_sequence):
        """
        Test valid extraction of coils when the structure has no secondary structures
        (i.e. no helices or sheets).
        """
        mock_structure.helices = []
        mock_structure.sheets = []

        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "A", 0, "SUBSEQ", "SUBSEQ", 1, 10, 6), 
            ("1A00", 2, "B", 0, "SUBSEQ", "SUBSEQ", 1, 8, 6)
        ]  

        # Expected output: 2 coils extracted
        assert len(result) == 2  
        # Expected output: entire sequence of chains A and B are extracted as coils 
        assert result == expected


    def test_insert_into_coil_table_unconfirmed_chain(self, mock_structure, mock_polymer_sequence, mock_chain_a, mock_chain_b):
        """
        Test valid extraction of coils when the whole chain is experimentally unconfirmed 
        (i.e. chain_object is None).
        """
        # chain_object is None
        mock_structure[0].find_chain.side_effect = lambda x: None 

        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=5)]]
        mock_chain_b.__getitem__.side_effect = [[MagicMock(label_seq=6)], [MagicMock(label_seq=7)]]  

        # mock start and end chains -- here, the helix and sheet are each in one chain
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
        ]

        # annotated subsequence is empty
        mock_polymer_sequence.get_chain_annotated_subsequence.return_value = ""
        # sequence contains unconfirmed residues 
        mock_polymer_sequence.contains_unconfirmed_residues.return_value = 1

        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "A", 1, "SUBSEQ", "", 6, 10, 6), 
            ("1A00", 2, "B", 1, "SUBSEQ", "", 1, 5, 6),
            ("1A00", 3, "B", 1, "SUBSEQ", "", 8, 8, 6)
        ]

        # Expected output: 3 coils (before and after helix/sheet)
        assert len(result) == 3
        assert result == expected

    # coil start > coil_end 
    def test_insert_into_coil_table_helix_at_chain_end(self, mock_structure, mock_chain_a, mock_chain_b, mock_polymer_sequence):
        """
        Test valid extraction of coils when a helix ends at the last residue of a chain. 
        """
        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=10)]]
        mock_chain_b.__getitem__.side_effect = [[MagicMock(label_seq=6)], [MagicMock(label_seq=7)]]                               

        # mock start and end chains -- here, the helix and sheet are each in one chain
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
        ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "B", 0, "SUBSEQ", "SUBSEQ", 1, 5, 6), 
            ("1A00", 2, "B", 0, "SUBSEQ", "SUBSEQ", 8, 8, 6)
        ]

        # Expected output: 2 coils (before and after helix/sheet)
        assert len(result) == 2
        assert result == expected
    

    def test_insert_into_coil_table_strand_at_chain_end(self, mock_structure, mock_chain_a, mock_chain_b, mock_polymer_sequence):
        """
        Test valid extraction of coils when a strand ends at the last residue of a chain. 
        """
        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=5)]]
        mock_chain_b.__getitem__.side_effect = [[MagicMock(label_seq=6)], [MagicMock(label_seq=8)]]                               

        # mock start and end chains -- here, the helix and sheet are each in one chain
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
        ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "A", 0, "SUBSEQ", "SUBSEQ", 6, 10, 6), 
            ("1A00", 2, "B", 0, "SUBSEQ", "SUBSEQ", 1, 5, 6)
        ]

        # Expected output: 2 coils (before and after helix/sheet)
        assert len(result) == 2
        assert result == expected

    
    def test_insert_into_coil_table_multiple_secondary_structures(self, mock_structure, mock_polymer_sequence, mock_chain_a, mock_chain_b):
        """
        Test valid extraction of coils with helices and sheets when a single chain 
        contains multiple secondary structures. 
        """
        # add another helix to the mock gemmi structure 
        helix2 = MagicMock()
        helix2.start.res_id.seqid.num = 7
        helix2.end.res_id.seqid.num = 7
        mock_structure.helices.append(helix2)

        # mock start and end ids of helices and strands 
        mock_chain_a.__getitem__.side_effect = [[MagicMock(label_seq=1)], [MagicMock(label_seq=5)],
                                                [MagicMock(label_seq=7)], [MagicMock(label_seq=7)]]
        mock_chain_b.__getitem__.side_effect = [[MagicMock(label_seq=6)], [MagicMock(label_seq=7)]]  

        # mock start and end chains -- here, two helices are in chain A, while the strand is in chain B.
        mock_structure[0].find_cra.side_effect = [
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_a), 
            MagicMock(chain=mock_chain_b), 
            MagicMock(chain=mock_chain_b)
        ]
        
        result = extract.insert_into_coil_table(mock_structure, MagicMock(), mock_polymer_sequence)
        expected = [
            ("1A00", 1, "A", 0, "SUBSEQ", "SUBSEQ", 6, 6, 6), 
            ("1A00", 2, "A", 0, "SUBSEQ", "SUBSEQ", 8, 10, 6), 
            ("1A00", 3, "B", 0, "SUBSEQ", "SUBSEQ", 1, 5, 6),
            ("1A00", 4, "B", 0, "SUBSEQ", "SUBSEQ", 8, 8, 6)
        ]

        # Expected output: 4 coils (before and after helix/sheet)
        assert len(result) == 4
        assert result == expected



