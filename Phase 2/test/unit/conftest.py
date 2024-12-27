import pytest
from unittest.mock import MagicMock
import gemmi 
from gemmi import cif, EntityType, PolymerType, EntityList
from sqlite3 import Cursor 

from polymer_sequence import PolymerSequence, Monomer
from extract import ComplexType
from table import Table
from attributes import Attributes

@pytest.fixture
def mock_structure():
    mock_structure = MagicMock(spec=gemmi.Structure)
    mock_structure.info = {"_entry.id": "1A00", "_struct.title": "mock_title", "_cell.Z_PDB": "1"}
    mock_structure.name = "mock_name"

    cell = MagicMock(spec=gemmi.UnitCell)
    cell.a, cell.b, cell.c = (1.0, 1.0, 1.0)
    cell.alpha, cell.beta, cell.gamma = (90.0, 90.0, 90.0)
    mock_structure.cell = cell

    mock_structure.spacegroup_hm = "P 1"

    return mock_structure

@pytest.fixture
def mock_doc():
    mock_doc = MagicMock(spec=cif.Document)
    block = MagicMock()
    mock_doc.sole_block.return_value = block

    # Mocking the loops in the CIF block
    block.find_loop.side_effect = lambda x: {
        "_pdbx_poly_seq_scheme.pdb_strand_id": ["A", "A", "B", "B", "B"],
        "_pdbx_poly_seq_scheme.entity_id": ["1", "1", "1", "1", "1"],
        "_pdbx_poly_seq_scheme.seq_id": ["1", "2", "3", "3"],
        "_pdbx_poly_seq_scheme.mon_id": ["ALA", "ARG", "ASN", "ASN"],
        "_pdbx_poly_seq_scheme.pdb_mon_id": ["ALA", "?", "ASN", "ASN"],
        "_pdbx_poly_seq_scheme.hetero": ["n", "n", "y", "y"]
    }[x]

    return mock_doc


@pytest.fixture
def fake_sequence_3to1():
    """Fixture to mock the sequence_3to1 function."""
    def sequence_3to1(monomer_sequence):
        mapping = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D"}
        return "".join(mapping[monomer] for monomer in monomer_sequence)
    return sequence_3to1

 
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
    mock_subchain.make_one_letter_sequence.return_value = "ARNDCQEGHIX"
    mock_subchain.length.return_value = 11

    return mock_subchain


@pytest.fixture
def mock_polymer():
    mock_polymer = MagicMock(spec=gemmi.ResidueSpan)

    # set up return value for len(mock_polymer)
    mock_polymer.__len__.return_value = 11 
    mock_polymer.make_one_letter_sequence.return_value = "ARNDCQEGHIX"
    mock_polymer.length.return_value = 11

    return mock_polymer


@pytest.fixture 
def mock_empty_polymer():
    mock_polymer = MagicMock(spec=gemmi.ResidueSpan)

    # set up return value for len(mock_polymer)
    mock_polymer.__len__.return_value = 0
    mock_polymer.make_one_letter_sequence.return_value = ""
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
def test_polymer_sequence():
    mock_doc = MagicMock(spec=cif.Document)
    test_polymer_sequence = PolymerSequence(mock_doc)
    
    test_polymer_sequence.one_letter_code = "ARNDCQEGHIX"

    test_polymer_sequence.chain_start_indices = {'A': 0, 'B': 0}
    test_polymer_sequence.chain_end_indices = {'A': 10}
    
    sequence = [Monomer("A", 1, 1, "ALA", "ALA", "n"), Monomer("A", 1, 2, "ARG", "ARG", "n"),
                Monomer("A", 1, 3, "ASN", "ASN", "n"), Monomer("A", 1, 4, "ASP", "ASP", "n"),
                Monomer("A", 1, 5, "CYS", "CYS", "n"), Monomer("A", 1, 6, "GLN", "GLN", "n"), 
                Monomer("A", 1, 7, "GLU", "GLU", "n"), Monomer("A", 1, 8, "GLY", "GLY", "n"), 
                Monomer("A", 1, 9, "HIS", "HIS", "n"), Monomer("A", 1, 10, "ILE", "ILE", "n"), 
                Monomer("A", 1, 11, "UNK", "UNK", "n")]
    test_polymer_sequence.sequence = sequence

    return test_polymer_sequence


@pytest.fixture
def mock_span():
    residues = []
    for label_seq in range(1, 11):  # Generate residues with label_seq 1 to 10
        mock_residue = MagicMock(spec=gemmi.Residue)
        mock_residue.label_seq = label_seq
        residues.append(mock_residue)
    return residues


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
def mock_table_schemas(mock_table, mock_coil_table):
    mock_table_schemas = [mock_table, mock_coil_table]
    return mock_table_schemas


@pytest.fixture
def mock_table():
    test_data = ("1A00", "data1", "data2")
    test_statement = "INSERT INTO main VALUES(?, ?, ?)"

    mock_table = MagicMock(spec=Table)
    mock_table.name = "main"
    mock_table.extract_data.return_value = [test_data]
    mock_table.insert_row.return_value = test_statement 

    return mock_table


@pytest.fixture
def mock_coil_table():
    test_data_1 = ("1A00", "data1", "data2")
    test_data_2 = ("1A00", "data3", "data4")
    test_statement = "INSERT INTO coils VALUES(?, ?, ?)"

    mock_table = MagicMock(spec=Table)
    mock_table.name = "coils"
    mock_table.extract_data.return_value = [test_data_1, test_data_2]
    mock_table.insert_row.return_value = test_statement 

    return mock_table


@pytest.fixture
def test_table():
    mock_attributes = MagicMock()
    mock_attributes.attribute_names = ("id", "a")
    mock_attributes.__str__.return_value = "(id VARCHAR, a FLOAT, PRIMARY KEY (id, a), FOREIGN KEY (id) REFERENCES main (id))"
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