import gemmi
from gemmi import cif, EntityType, PolymerType
import sqlite3
import attr
from attr import ComplexType

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

def has_nonpolymer(struct: gemmi.Structure):
    for entity in struct.entities:
        if entity.entity_type == EntityType.NonPolymer:
            return True
    return False

def has_antibody(struct: gemmi.Structure):
    return "antibody" in struct.name

def get_complex_type(struct: gemmi.Structure):
    # We use bitwise operations to quickly get the complex type

    # these values are 0 if compound doesn't have
    # a saccharide/nucleic acid/peptide, and nonzero if it does
    has_saccharide = has_nucleic_acid = has_peptide = 0

    num_peptides = 0
    peptide_entity = None

    for entity in struct.entities:
        if entity.entity_type == EntityType.Unknown:
            return ComplexType.Other
        elif entity.entity_type == EntityType.Branched:
            has_saccharide = 0b100
        elif entity.entity_type == EntityType.Polymer:
            if entity.polymer_type in [PolymerType.Other, PolymerType.Unknown, PolymerType.CyclicPseudoPeptide]:
                return ComplexType.Other
            elif entity.polymer_type in [PolymerType.Dna, PolymerType.Rna, PolymerType.DnaRnaHybrid, PolymerType.Pna]:
                has_nucleic_acid = 0b010
            elif entity.polymer_type in [PolymerType.SaccharideD, PolymerType.SaccharideL]:
                has_saccharide = 0b100
            else: # entity.polymer_type in [PolymerType.PeptideL, PolymerType.PeptideD]
                num_peptides += 1
                has_peptide = 0b001
                peptide_entity = entity
        
        # regular integer addition gives the same result, but bitwise or makes the intention clearer
        pending_complex_type = ComplexType(has_peptide | has_nucleic_acid | has_saccharide)

        # Check if compound is actually a complex protein or proteinmer
        if pending_complex_type == ComplexType.SingleProtein: # Compound has polypeptides, and no other polymers
            if num_peptides > 1:
                return ComplexType.ComplexProtein
            else: # num_peptides == 1
                if len(peptide_entity.subchains) > 1:
                    return ComplexType.Proteinmer
        
        return pending_complex_type

def insert_into_main_table(struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    complex_type = get_complex_type(struct)
    chains = [chain.name for chain in struct[0]]
    cell = struct.cell
    z_value = None
    if "_cell.Z_PDB" in struct.info:
        z_value = struct.info["_cell.Z_PDB"]
    spacegroup = struct.spacegroup_hm
    cur.execute("INSERT INTO " + attr.main_table + " VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (id, complex_type.name, has_nonpolymer(struct), has_antibody(struct), ' '.join(chains),
                 spacegroup, z_value, cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma))
        
def insert_into_entity_table(doc: cif.Document, struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    names = block.find_loop("_entity.pdbx_description")
    if len(names) == 0:
        name = block.find_value("_entity.pdbx_description")
        names = [name]
    for index, entity in enumerate(struct.entities):
        cur.execute("INSERT INTO " + attr.entity_table + " VALUES(?, ?, ?, ?, ?, ?)",
                    (id, entity.name, names[index], entity.entity_type.name,
                     entity.polymer_type.name, ' '.join(entity.subchains)))
        
def insert_into_subchain_table(struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for entity in struct.entities:
        if entity.polymer_type in [PolymerType.PeptideD, PolymerType.PeptideL]:
            for subchain_name in entity.subchains:
                subchain = struct[0].get_subchain(subchain_name)
                parent_chain = struct[0].get_parent_of(subchain[0]).name
                cur.execute("INSERT INTO " + attr.subchain_table + " VALUES(?, ?, ?, ?, ?, ?)",
                            (id, entity.name, subchain.subchain_id(), parent_chain,
                             subchain.make_one_letter_sequence(), subchain.length()))
        
def insert_into_helix_table(struct: gemmi.Structure, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for helix in struct.helices:
        chain = struct[0].find_cra(helix.start).chain
        sequence = sequence_3to1(chain.whole().extract_sequence(),
                                 helix.start.res_id.seqid.num - 1, helix.end.res_id.seqid.num)
        cur.execute("INSERT INTO " + attr.helix_table + " VALUES(?, ?, ?, ?)",
                    (id, chain.name, sequence, helix.length))
    
        
def insert_into_all_tables(path: str, cur: sqlite3.Cursor):
    #try:
        struct = gemmi.read_structure(path)
        doc = cif.read_file(path)
        
        insert_into_main_table(struct, cur)
        insert_into_entity_table(doc, struct, cur)
        insert_into_subchain_table(struct, cur)
        insert_into_helix_table(struct, cur)
    #except Exception as error:
    #    print("Error at " + path)
    #    print(error)