import gemmi
from gemmi import cif, EntityType, PolymerType
import sqlite3
import attr
from attr import ComplexType

# Extract a substring from the one-letter encoding of a given sequence of residues.
# Start and end indices should be given in the one-indexed sequence address form, so
# the first residue in the sequence has index 1, and the function returns the
# substring from the start index to end index inclusive.
# If the start index is greater than the end index, then the output sequence is
# given in reverse of its input sequence.
#
# This function is intended to work for peptide polymer residue spans only, and is highly
# susceptible to errors for other types of spans.
def one_letter_sequence(chain: gemmi.Chain, start: int = -1, end: int = -1) -> str:
    sequence = chain.whole().make_one_letter_sequence()
    first_auth_index = chain[0].seqid
    first_index = chain[str(first_auth_index)].auth_seq_id_to_label(first_auth_index)
    last_auth_index = chain[-1].seqid
    last_index = chain[str(last_auth_index)].auth_seq_id_to_label(last_auth_index)
    # if no values given for start and end, assign them to start and end of span
    if start == -1:
        start = first_index
    if end == -1:
        end = last_index
    start -= first_index
    end -= first_index

    if end >= start:
        return sequence[start:end + 1]
    if end == 0:
        return sequence[start::-1]
    return sequence[start:end-1:-1]

# For the beta sheets table. Produces a sequence of the series of
# parallel and antiparallel bonds between strands
def sense_sequence(sheet: gemmi.Sheet) -> str:
    encoding = {1: 'P', -1: 'A', 0: ''}
    return ''.join([encoding[strand.sense] for strand in sheet.strands])

# Determines the type of complex of structure based on the types of its entities
def get_complex_type(struct: gemmi.Structure) -> ComplexType:
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

# Insert given data into the given table
def insert_into_table(cur: sqlite3.Cursor, table_name: str, data):
    args = ', '.join(['?' for i in range(len(data))])
    query = f'INSERT INTO {table_name} VALUES({args})'
    cur.execute(query, data)

def insert_into_main_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    source_org = block.find_value("_entity_src_gen.pdbx_gene_src_scientific_name")
    complex_type = get_complex_type(struct)
    chains = [chain.name for chain in struct[0]]
    cell = struct.cell
    z_value = None
    if "_cell.Z_PDB" in struct.info:
        z_value = struct.info["_cell.Z_PDB"]
    spacegroup = struct.spacegroup_hm
    data = (id, complex_type.name, source_org, ' '.join(chains), spacegroup, z_value,
                 cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    insert_into_table(cur, attr.main_table[0], data)
        
def insert_into_entity_table(struct: gemmi.Structure, doc: cif.Document, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    names = block.find_loop("_entity.pdbx_description")
    if len(names) == 0:
        name = block.find_value("_entity.pdbx_description")
        names = [name]
    for index, entity in enumerate(struct.entities):
        data = (id, entity.name, names[index], entity.entity_type.name,
                     entity.polymer_type.name, ' '.join(entity.subchains))
        insert_into_table(cur, attr.entity_table[0], data)
        
def insert_into_subchain_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for entity in struct.entities:
        if entity.polymer_type in [PolymerType.PeptideD, PolymerType.PeptideL]:
            for subchain_name in entity.subchains:
                subchain = struct[0].get_subchain(subchain_name)
                if len(subchain) == 0:
                    continue
                parent_chain = struct[0].get_parent_of(subchain[0])
                start_auth_pos = subchain[0].seqid
                end_auth_pos = subchain[-1].seqid
                start_pos = parent_chain.whole()[str(start_auth_pos)].auth_seq_id_to_label(start_auth_pos)
                end_pos = parent_chain.whole()[str(end_auth_pos)].auth_seq_id_to_label(end_auth_pos)
                data = (id, entity.name, subchain.subchain_id(), parent_chain.name,
                        subchain.make_one_letter_sequence(), start_pos, end_pos, subchain.length())
                insert_into_table(cur, attr.subchain_table[0], data)

def insert_into_chain_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for chain in struct[0]:
        start_auth_pos = chain.whole()[0].seqid
        end_auth_pos = chain.whole()[-1].seqid
        start_pos = chain.whole().auth_seq_id_to_label(start_auth_pos)
        end_pos = chain.whole().auth_seq_id_to_label(end_auth_pos)
        data = (id, chain.name, ' '.join([subchain.subchain_id() for subchain in chain.subchains()]),
                chain.whole().make_one_letter_sequence(), start_pos, end_pos, chain.whole().length())
        insert_into_table(cur, attr.chain_table[0], data)
        
def insert_into_helix_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for helix in struct.helices:
        chain = struct[0].find_cra(helix.start).chain
        end_chain = struct[0].find_cra(helix.end).chain
        start_auth_pos = helix.start.res_id.seqid
        end_auth_pos = helix.end.res_id.seqid

        if (chain != end_chain):
            start_pos = chain.whole().auth_seq_id_to_label(start_auth_pos)
            end_pos = end_chain.whole().auth_seq_id_to_label(end_auth_pos)
            length = helix.length
            sequence = "MULTIPLE CHAINS ERROR"
            data = (id, chain.name + ' ' + end_chain.name, sequence, start_pos, end_pos, length)
        else:
            start_pos = chain.whole().auth_seq_id_to_label(start_auth_pos)
            end_pos = chain.whole().auth_seq_id_to_label(end_auth_pos)
            direction = 1 if start_pos <= end_pos else -1
            length = end_pos - start_pos + direction
            sequence = one_letter_sequence(chain, start_pos, end_pos)
            data = (id, chain.name, sequence, start_pos, end_pos, length)

        insert_into_table(cur, attr.helix_table[0], data)
        
def insert_into_sheet_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for sheet in struct.sheets:
        data = (id, sheet.name, len(sheet.strands), sense_sequence(sheet))
        insert_into_table(cur, attr.sheet_table[0], data)

def insert_into_strand_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for sheet in struct.sheets:
        for strand in sheet.strands:
            chain = struct[0].find_cra(strand.start).chain
            start_auth_pos = strand.start.res_id.seqid
            end_auth_pos = strand.end.res_id.seqid
            start_pos = chain.whole().auth_seq_id_to_label(start_auth_pos)
            end_pos = chain.whole().auth_seq_id_to_label(end_auth_pos)
            direction = 1 if start_pos <= end_pos else -1
            length = end_pos - start_pos + direction
            sequence = one_letter_sequence(chain, start_pos, end_pos)
            data = (id, sheet.name, strand.name, chain.name, sequence, start_pos, end_pos, length)
            insert_into_table(cur, attr.strand_table[0], data)

# Inserts the data of a given file into all tables
def insert_into_all_tables(path: str, cur: sqlite3.Cursor):

    insert_functions = [insert_into_main_table, insert_into_entity_table, insert_into_subchain_table,
                        insert_into_chain_table, insert_into_helix_table, insert_into_sheet_table,
                        insert_into_strand_table]
    
    try:
        struct = gemmi.read_structure(path)
        doc = cif.read(path)
        res = cur.execute("SELECT entry_id FROM " + attr.main_table[0]\
                            + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
        if not res.fetchone(): # if there is no row in the main table with such entry ID
            for function in insert_functions:
                function(struct, doc, cur)
    except Exception as error:
        print("Error at " + path)
        print(error)