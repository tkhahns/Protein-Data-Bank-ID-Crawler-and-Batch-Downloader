import gemmi
from gemmi import cif
import sqlite3
import attr

three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}

def letter_code_3to1(polymer: str) -> str:
    if (polymer in three_to_one.keys()):
        return three_to_one[polymer]
    return '?'

def sequence_3to1(sequence: list[str], start: int = 0, end: int = -1) -> str:
    return ''.join([letter_code_3to1(polymer) for polymer in sequence[start:end]])

def insert_into_main_table(doc: cif.Document, struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    proteins = block.find_value("_refine_hist.pdbx_number_atoms_protein")
    acids = block.find_value("_refine_hist.pdbx_number_atoms_nucleic_acid")
    cell = struct.cell
    z_value = None
    if "_cell.Z_PDB" in struct.info:
        z_value = struct.info["_cell.Z_PDB"]
    spacegroup = struct.spacegroup_hm
    cur.execute("INSERT INTO " + attr.main_table + " VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (id, proteins, acids, spacegroup, z_value,
                cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma))
        
def insert_into_entity_table(doc: cif.Document, struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    names = block.find_loop("_entity.pdbx_description")
    if names is None:
        name = block.find_value("_entity.pdbx_description")
        names = [name]
    for index, entity in enumerate(struct.entities):
        cur.execute("INSERT INTO " + attr.entity_table + " VALUES(?, ?, ?, ?, ?)",
                    (id, names[index], str(entity.entity_type),
                    ' '.join(entity.subchains), sequence_3to1(entity.full_sequence)))
        
def insert_into_chain_table(struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for chain in struct[0]:
        sequence = sequence_3to1(chain.whole().extract_sequence())
        cur.execute("INSERT INTO " + attr.chain_table + " VALUES(?, ?, ?)",
                    (id, chain.name, sequence))
        
def insert_into_helix_table(struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for helix in struct.helices:
        chain = struct[0].find_cra(helix.start).chain
        sequence = sequence_3to1(chain.whole().extract_sequence(),
                                 helix.start.res_id.seqid.num - 1, helix.end.res_id.seqid.num)
        cur.execute("INSERT INTO " + attr.helix_table + " VALUES(?, ?, ?)",
                    (id, chain.name, sequence))
        
def insert_into_all_tables(path: str, cur: sqlite3.Cursor):
    struct = gemmi.read_structure(path)
    doc = cif.read_file(path)
    block = doc.sole_block()
    proteins = block.find_value("_refine_hist.pdbx_number_atoms_protein")
    if proteins is not None and proteins == 0:
        return
    
    insert_into_main_table(doc, struct, cur)
    insert_into_entity_table(doc, struct, cur)
    insert_into_chain_table(struct, cur)
    insert_into_helix_table(struct, cur)