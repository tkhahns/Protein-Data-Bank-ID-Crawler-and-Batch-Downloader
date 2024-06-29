__author__ = "Siwan Li"
"""
This script contains all the functions needed for extracting information out of mmCIF files.
Important things to note:
- ResidueSpan.length() is often used over len(ResidueSpan), as the former counts the length of
  the residue span (so only one residue  out of a set of microheterogeneities are counted),
  whereas the latter counts how many residue elements are contained within the span
  (so experimentally unconfirmed residues are excluded, and all residues of a set of
  microheterogeneities are counted). Note that experimentally unconfirmed residues are
  excluded from the length.
- While either the primary sequence ID or the author sequence ID could be used, the primary
  sequence ID is consistently used as it is always strictly increasing along the span. The author
  sequence ID may decrease at certain points, as many residues may have the same numerical
  author sequence ID, but have different 'icode's (see gemmi.SeqId.icode).
"""

import gemmi
from gemmi import cif, EntityType, PolymerType
import sqlite3
import attr
from attr import ComplexType

def one_letter_sequence(span: gemmi.ResidueSpan, start: int = -1, end: int = -1) -> str:
    """
    Extract a substring from the one-letter encoding of a given sequence of residues.
    Start and end indices should be given in label sequence ID form,
    i.e. the integer given from Residue.label_seq, and the function returns the
    substring from the start index to end index inclusive.
    If the start index is greater than the end index, then the output sequence is
    given in reverse of its input sequence.
    """
    sequence = span.make_one_letter_sequence()
    first_label = span[0].label_seq
    last_label = span[-1].label_seq
    # if no values given for start and end, assign them to start and end of span
    if start == -1:
        start = first_label
    start, string_start = get_label_and_string_index(span, start)
    if end == -1:
        end = last_label
    end, string_end = get_label_and_string_index(span, end)

    if end >= start:
        return sequence[string_start:string_end + 1], end - start + 1
    if end == 0:
        return sequence[string_start::-1], -start - 1
    return sequence[string_start:string_end-1:-1], end - start - 1

def get_label_and_string_index(span: gemmi.ResidueSpan, target_label: int) -> int:
    """
    Given some residue span and a target label ID, this returns the 'label index'
    (i.e. the index in the span that calls the residue with the desired label sequence ID),
    and the 'string index' (i.e. the index in the span's one letter sequence that
    calls the character of the corresponding residue with the desired label sequence ID).
    """
    first_conformer_span = list(span.first_conformer())
    target_index = binary_search(first_conformer_span, 0, len(first_conformer_span) - 1, target_label)
    sequence = span.make_one_letter_sequence()
    dashes = sequence[:target_index + 1].count('-')
    string_index = target_index
    while (dashes != 0):
        next_index = string_index + dashes
        dashes = sequence[string_index + 1:next_index + 1].count('-')
        string_index = next_index
    return target_index, string_index

def binary_search(span: gemmi.ResidueSpan, left_index: int, right_index: int, target_label: int) -> int:
    """
    Helper method for finding the index of a desired element in a residue span by means of binary search.
    """
    if right_index < left_index:
        raise Exception("Couldn't find index")
    centre_index = (left_index + right_index) // 2
    centre_label = span[centre_index].label_seq
    if centre_label == target_label:
        return centre_index
    difference = target_label - centre_label
    next_index = centre_index + difference
    # ensure that next_index are between the left_index and right_index
    next_index = max(left_index, next_index)
    next_index = min(right_index, next_index)
    next_label = span[next_index].label_seq
    if next_label == target_label:
        return next_index
    
    if target_label < centre_label:
        return binary_search(span, left_index, centre_index - 1, target_label)
    # centre_label < target_label
    return binary_search(span, centre_index + 1, right_index, target_label)
        
def sense_sequence(sheet: gemmi.Sheet) -> str:
    """
    For the beta sheets table. Produces a sequence of the series of
    parallel and antiparallel bonds between strands.
    """
    encoding = {1: 'P', -1: 'A', 0: ''}
    return ''.join([encoding[strand.sense] for strand in sheet.strands])

def get_complex_type(struct: gemmi.Structure) -> ComplexType:
    """
    Determines the type of complex of the protein structure based on the types of its entities.
    """
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

def insert_into_table(cur: sqlite3.Cursor, table_name: str, data):
    """
    Inserts the given data into the given table
    """
    args = ', '.join(['?' for i in range(len(data))])
    query = f'INSERT INTO {table_name} VALUES({args})'
    cur.execute(query, data)

def insert_into_main_table(struct: gemmi.Structure, doc: cif.Document, cur: sqlite3.Cursor):
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
                start_pos = subchain[0].label_seq
                end_pos = subchain[-1].label_seq
                data = (id, entity.name, subchain.subchain_id(), parent_chain.name,
                        subchain.make_one_letter_sequence(), start_pos, end_pos, subchain.length())
                insert_into_table(cur, attr.subchain_table[0], data)

def insert_into_chain_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for chain in struct[0]:
        if len(chain.get_polymer()) == 0:
            start_pos = end_pos = None
        else:
            start_pos = chain.get_polymer()[0].label_seq
            end_pos = chain.get_polymer()[-1].label_seq
        sequence = chain.get_polymer().make_one_letter_sequence()
        data = (id, chain.name, ' '.join([subchain.subchain_id() for subchain in chain.subchains()]),
                sequence, start_pos, end_pos, chain.get_polymer().length())
        insert_into_table(cur, attr.chain_table[0], data)
        
def insert_into_helix_table(struct: gemmi.Structure, doc, cur: sqlite3.Cursor):
    id = struct.info["_entry.id"]
    for helix in struct.helices:
        chain = struct[0].find_cra(helix.start).chain
        end_chain = struct[0].find_cra(helix.end).chain
        start_auth_label = str(helix.start.res_id.seqid.num) + helix.start.res_id.seqid.icode
        end_auth_label = str(helix.end.res_id.seqid.num) + helix.end.res_id.seqid.icode
        start_pos = chain[start_auth_label][0].label_seq
        end_pos = end_chain[end_auth_label][0].label_seq

        if (chain != end_chain):
            length = helix.length
            sequence = "MULTIPLE CHAINS ERROR"
            data = (id, chain.name + ' ' + end_chain.name, sequence, start_pos, end_pos, length)
        else:
            sequence, length = one_letter_sequence(chain.get_polymer(), start_pos, end_pos)
            data = (id, chain.name, sequence, start_pos, end_pos, helix.length)

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
            end_chain = struct[0].find_cra(strand.end).chain
            start_auth_label = str(strand.start.res_id.seqid.num) + strand.start.res_id.seqid.icode
            end_auth_label = str(strand.end.res_id.seqid.num) + strand.end.res_id.seqid.icode
            start_pos = chain[start_auth_label][0].label_seq
            end_pos = end_chain[end_auth_label][0].label_seq

            if (chain != end_chain):
                sequence = "MULTIPLE CHAINS ERROR"
                data = (id, sheet.name, strand.name, chain.name + ' ' + end_chain.name, sequence, start_pos, end_pos, -1)
            sequence, length = one_letter_sequence(chain.get_polymer(), start_pos, end_pos)
            data = (id, sheet.name, strand.name, chain.name, sequence, start_pos, end_pos, length)
            insert_into_table(cur, attr.strand_table[0], data)

def insert_into_all_tables(path: str, cur: sqlite3.Cursor):
    """
    Inserts the data of a given file into all tables
    """

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
        print(struct.name)
        print(error)