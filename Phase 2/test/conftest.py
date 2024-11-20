import pytest
from unittest.mock import MagicMock
import gemmi 
from gemmi import cif, EntityType, PolymerType, EntityList
from sqlite3 import Cursor 

from polymer_sequence import PolymerSequence
from extract import ComplexType
from table import Table
from attributes import Attributes

@pytest.fixture
def mock_structure():
    mock_structure = MagicMock(spec=gemmi.Structure)
    mock_structure.info = {'_entry.id': '1A00', '_struct.title': 'mock_title', '_cell.Z_PDB': '1'}
    mock_structure.name = "mock_name"

    cell = MagicMock(spec=gemmi.UnitCell)
    cell.a, cell.b, cell.c = (1.0, 1.0, 1.0)
    cell.alpha, cell.beta, cell.gamma = (90.0, 90.0, 90.0)
    mock_structure.cell = cell

    mock_structure.spacegroup_hm = 'P 1'

    return mock_structure

@pytest.fixture
def mock_doc():
    mock_doc = MagicMock(spec=cif.Document)
    return mock_doc
 
@pytest.fixture
def mock_entity():
    def create_mock_entity(entity_type, polymer_type = None, subchains = []):
        mock_entity = gemmi.Entity("")
        mock_entity.name = '1'
        mock_entity.entity_type = entity_type
        if polymer_type:
            mock_entity.polymer_type = polymer_type
        mock_entity.subchains = subchains

        return mock_entity
    
    return create_mock_entity


@pytest.fixture
def mock_entities(mock_entity):
    def create_entities_with_complex_type(complex_type):
        entities = EntityList()

        # arguments for creating different entities given a complex type
        entity_args = {
            ComplexType.Other: [
                (EntityType.Unknown,)
            ],
            ComplexType.NucleicAcid : [
                (EntityType.Polymer, PolymerType.Dna), 
                (EntityType.Polymer, PolymerType.Rna), 
                (EntityType.Polymer, PolymerType.DnaRnaHybrid),
                (EntityType.Polymer, PolymerType.Pna)
            ],
            ComplexType.Saccharide: [
                (EntityType.Branched,),
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.SaccharideL)
            ],
            ComplexType.SingleProtein: [
                (EntityType.Polymer, PolymerType.PeptideL, ["A"])
            ],
            ComplexType.Proteinmer: [
                (EntityType.Polymer, PolymerType.PeptideL, ["A", "B"])
            ], 
            ComplexType.ComplexProtein: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.PeptideL)
            ],
            ComplexType.ProteinNA: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.Dna)
            ],
            ComplexType.ProteinSaccharide: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.SaccharideD)
            ],
            ComplexType.SaccharideNA: [
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.Dna)
            ], 
            ComplexType.ProteinSaccharideNA: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.Dna)
            ]    
        }

        for args in entity_args[complex_type]:
            entities.append(mock_entity(*args))

        return entities 

    return create_entities_with_complex_type


@pytest.fixture
def mock_subchain():
    mock_subchain = MagicMock(spec=gemmi.ResidueSpan)
    mock_subchain.subchain_id.return_value = 'A'
    mock_subchain.make_one_letter_sequence.return_value = 'ARNDCQEGHIX'
    mock_subchain.length.return_value = 11

    return mock_subchain


@pytest.fixture
def mock_polymer():
    mock_polymer = MagicMock(spec=gemmi.ResidueSpan)

    # set up return value for len(mock_polymer)
    mock_polymer.__len__.return_value = 11 
    mock_polymer.make_one_letter_sequence.return_value = 'ARNDCQEGHIX'
    mock_polymer.length.return_value = 11

    return mock_polymer


@pytest.fixture 
def mock_empty_polymer():
    mock_polymer = MagicMock(spec=gemmi.ResidueSpan)

    # set up return value for len(mock_polymer)
    mock_polymer.__len__.return_value = 0
    mock_polymer.make_one_letter_sequence.return_value = ''  
    mock_polymer.length.return_value = 0  

    return mock_polymer


@pytest.fixture
def mock_chain(mock_subchain, mock_polymer):
    mock_chain = MagicMock(spec=gemmi.Chain)
    mock_chain.name = 'A'

    mock_chain.get_polymer.return_value = mock_polymer
    mock_chain.subchains.return_value = [mock_subchain, mock_subchain]

    return mock_chain


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
    mock_polymer_sequence = PolymerSequence(mock_doc)
    
    mock_polymer_sequence.one_letter_code = 'ARNDCQEGHIX'

    mock_polymer_sequence.chain_start_indices = {'A': 0, 'B': 0}
    mock_polymer_sequence.chain_end_indices = {'A': 10}

    sequence = [(0, 0, 1), (0, 0, 3), (0, 0, 5), (0, 0, 7), (0, 0, 9)]
    mock_polymer_sequence.sequence = sequence

    return mock_polymer_sequence


@pytest.fixture
def mock_helix():
    mock_helix = MagicMock(spec=gemmi.Helix)
    mock_helix.start = MagicMock(spec=gemmi.AtomAddress)
    mock_helix.end = MagicMock(spec=gemmi.AtomAddress)
    start = mock_helix.start
    end = mock_helix.end

    start.res_id.seqid.num = 1
    start.res_id.seqid.icode = ''
    end.res_id.seqid.num = 11
    end.res_id.seqid.icode = ''

    mock_helix.length = 11

    return mock_helix


@pytest.fixture 
def mock_sheet():
    mock_sheet = MagicMock(spec=gemmi.Sheet)
    mock_sheet.name = 'A'

    return mock_sheet


@pytest.fixture
def mock_strand():
    mock_strand = MagicMock(spec=gemmi.Sheet.Strand)
    mock_strand.name = '1'
    mock_strand.start = MagicMock(spec=gemmi.AtomAddress)
    mock_strand.end = MagicMock(spec=gemmi.AtomAddress)
    start = mock_strand.start
    end = mock_strand.end

    start.res_id.seqid.num = 1
    start.res_id.seqid.icode = ''
    end.res_id.seqid.num = 11
    end.res_id.seqid.icode = ''

    return mock_strand


@pytest.fixture
def mock_cursor():
    mock_cursor = MagicMock(spec=Cursor)
    return mock_cursor


@pytest.fixture
def mock_table():
    TEST_DATA = ('1A00', 'data1', 'data2')
    TEST_STATEMENT = "INSERT INTO main VALUES(?, ?, ?)"

    mock_table = MagicMock(spec=Table)
    mock_table.name = "main"
    mock_table.extract_data.return_value = [TEST_DATA]
    mock_table.insert_row.return_value = TEST_STATEMENT 

    return mock_table


@pytest.fixture
def test_table():
    mock_attributes = MagicMock()
    mock_attributes.__str__.return_value = "(id VARCHAR, a FLOAT,\
        PRIMARY KEY(id, a), FOREIGN KEY (id) REFERENCES\
                main(id))"
    mock_attributes.match_columns.return_value = "col1 = value1, col2 = value2"
    mock_attributes.match_primary_keys.return_value = "id = 1"

    mock_extractor = MagicMock(return_value=[("test_id", "data1")])

    return Table("test_table", mock_attributes, mock_extractor)


@pytest.fixture
def test_attributes():
    test_attributes_pairs = [("id", "VARCHAR"), ("a", "FLOAT")]
    test_primary_keys = ["id", "a"]
    test_foreign_keys = {"id": ("main", "id")}

    return Attributes(test_attributes_pairs, test_primary_keys, test_foreign_keys)