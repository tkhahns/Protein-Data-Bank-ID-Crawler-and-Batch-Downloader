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

from typing import NewType
import gemmi
from gemmi import cif, EntityType, PolymerType
from polymer_sequence import PolymerSequence
from enum import Enum

MainData = NewType("MainData", tuple[str, str, str, str, str, str, str, int, float, float, float, float, float, float])
ExperimentalData = NewType("ExperimentalData", tuple[float, float, str, str, str, str, float, float])
EntityData = NewType("EntityData", tuple[str, str, str, str, str, str])
ChainData = NewType("ChainData", tuple[str, str, str, int, str, str, int, int, int, int])
SubchainData = NewType("SubchainData", tuple[str, str, str, int, str, str, int, int, int])
HelixData = NewType("HelixData", tuple[str, str, int, str, int, int, int])
SheetData = NewType("SheetData", tuple[str, str, int, str])
StrandData = NewType("StrandData", tuple[str, str, str, str, int, str, int, int, int])
CoilData = NewType("CoilData", tuple[str, int, str, int, str, str, int, int, int])

# Possible types of a complex, based on their entities
ComplexType = Enum('ComplexType', ['Other', 'SingleProtein',
                 'NucleicAcid', 'ProteinNA', 'Saccharide',
                 'ProteinSaccharide', 'SaccharideNA', 'ProteinSaccharideNA',
                 'Proteinmer', 'ComplexProtein'], start=0)
        
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

def insert_into_main_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> MainData:
    id = struct.info["_entry.id"]
    struct_title = struct.info["_struct.title"]
    block = doc.sole_block()
    source_org = block.find_value("_entity_src_gen.pdbx_gene_src_scientific_name")
    if source_org is None:
        source_org = ''
    source_org = source_org.strip("'")
    revision_date = block.find_value("_pdbx_audit_revision_history.revision_date")
    if revision_date is None:
        revision_date = block.find_loop("_pdbx_audit_revision_history.revision_date")[-1]
    complex_type = get_complex_type(struct)
    chains = [chain.name for chain in struct[0]]
    cell = struct.cell
    z_value = ''
    if "_cell.Z_PDB" in struct.info:
        z_value = int(struct.info["_cell.Z_PDB"])
    spacegroup = struct.spacegroup_hm
    return [(id, complex_type.name, struct_title, source_org, revision_date, ' '.join(chains),
             spacegroup, z_value, cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)]

def insert_into_experimental_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> ExperimentalData:
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    matthews_coefficient = block.find_value("_exptl_crystal.density_Matthews")
    percent_solvent_content = block.find_value("_exptl_crystal.density_percent_sol")
    crystal_growth_method = block.find_value("_exptl_crystal_grow.method")
    crystal_growth_proc = block.find_value("_exptl_crystal_grow.pdbx_details")
    crystal_growth_apparatus = block.find_value("_exptl_crystal_grow.apparatus")
    crystal_growth_atmosphere = block.find_value("_exptl_crystal_grow.atmosphere")
    crystal_growth_pH = block.find_value("_exptl_crystal_grow.pH")
    crystal_growth_temp = block.find_value("_exptl_crystal_grow.temp")
    data = [id, matthews_coefficient, percent_solvent_content, crystal_growth_method, crystal_growth_proc,\
            crystal_growth_apparatus, crystal_growth_atmosphere, crystal_growth_pH, crystal_growth_temp]
    for i in range(len(data)):
        if data[i] is None:
            data[i] = ''
    return [tuple(data)]
        
def insert_into_entity_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> EntityData:
    data = []
    id = struct.info["_entry.id"]
    block = doc.sole_block()
    names = block.find_loop("_entity.pdbx_description")
    if len(names) == 0:
        name = block.find_value("_entity.pdbx_description")
        names = [name]
    for index, entity in enumerate(struct.entities):
        if names[index] is None:
            names[index] = ''
        data.append((id, entity.name, names[index].strip("'"), entity.entity_type.name, entity.polymer_type.name,
                     ' '.join(entity.subchains)))
    return data
        
def insert_into_subchain_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> SubchainData:
    data = []
    id = struct.info["_entry.id"]
    for entity in struct.entities:
        if entity.polymer_type in [PolymerType.PeptideD, PolymerType.PeptideL]:
            for subchain_name in entity.subchains:
                subchain = struct[0].get_subchain(subchain_name)
                if len(subchain) == 0:
                    continue
                parent_chain = struct[0].get_parent_of(subchain[0])
                start_id = subchain[0].label_seq
                end_id = subchain[-1].label_seq
                annotated_sequence = subchain.make_one_letter_sequence()
                unannotated_sequence = sequence.get_chain_subsequence(parent_chain.name, start_id, end_id)[0]
                data.append((id, entity.name, subchain.subchain_id(), parent_chain.name,
                        unannotated_sequence, annotated_sequence, start_id, end_id, subchain.length()))
    return data

def insert_into_chain_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> ChainData:
    data = []
    id = struct.info["_entry.id"]
    for chain in struct[0]:
        if len(chain.get_polymer()) == 0:
            start_id = end_id = unconfirmed = None
        else:
            start_id = sequence.get_chain_start_id(chain.name)
            end_id = sequence.get_chain_end_id(chain.name)
            unconfirmed = sequence.contains_unconfirmed_residues(chain.name, start_id, end_id)
        subchains = ' '.join([subchain.subchain_id() for subchain in chain.subchains()])
        annotated_sequence = chain.get_polymer().make_one_letter_sequence()
        unannotated_sequence = sequence.get_chain_sequence(chain.name)
        author_start_id = chain.get_polymer()[0].seqid.num if len(chain.get_polymer()) > 0 else None
        author_end_id = chain.get_polymer()[-1].seqid.num if len(chain.get_polymer()) > 0 else None
        data.append((id, chain.name, subchains, unconfirmed,
                unannotated_sequence, annotated_sequence, start_id, end_id,
                chain.get_polymer().length(), author_start_id, author_end_id))
    return data
        
def insert_into_helix_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> HelixData:
    data = []
    id = struct.info["_entry.id"]
    for index, helix in enumerate(struct.helices):
        helix_sequence = sequence.get_helix_sequence(helix, struct)
        chain = struct[0].find_cra(helix.start).chain
        end_chain = struct[0].find_cra(helix.end).chain
        start_auth_label = str(helix.start.res_id.seqid.num) + helix.start.res_id.seqid.icode
        end_auth_label = str(helix.end.res_id.seqid.num) + helix.end.res_id.seqid.icode
        start_id = chain[start_auth_label][0].label_seq
        end_id = end_chain[end_auth_label][0].label_seq
        if chain != end_chain:
            chain_names = chain.name + ' ' + end_chain.name
            data.append((id, index + 1, chain_names, helix_sequence, start_id, end_id, helix.length))
        else:
            data.append((id, index + 1, chain.name, helix_sequence, start_id, end_id, helix.length))
    return data
        
def insert_into_sheet_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> SheetData:
    data = []
    id = struct.info["_entry.id"]
    for sheet in struct.sheets:
        data.append((id, sheet.name, len(sheet.strands), sense_sequence(sheet)))
    return data

def insert_into_strand_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> StrandData:
    data = []
    id = struct.info["_entry.id"]
    for sheet in struct.sheets:
        for strand in sheet.strands:
            strand_sequence, length = sequence.get_strand_sequence(strand, struct)
            chain = struct[0].find_cra(strand.start).chain
            end_chain = struct[0].find_cra(strand.end).chain
            start_auth_label = str(strand.start.res_id.seqid.num) + strand.start.res_id.seqid.icode
            end_auth_label = str(strand.end.res_id.seqid.num) + strand.end.res_id.seqid.icode
            start_id = chain[start_auth_label][0].label_seq
            end_id = end_chain[end_auth_label][0].label_seq
            data.append((id, sheet.name, strand.name, chain.name,
                         strand_sequence, start_id, end_id, length))
    return data

def insert_into_coil_table(struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence) -> CoilData:
    data = []
    secondary_structures = []
    id = struct.info["_entry.id"]

    # We first gather all the helix and sheet secondary structures to 'remove' from the polymer sequences.
    for index, helix in enumerate(struct.helices):
        chain = struct[0].find_cra(helix.start).chain
        end_chain = struct[0].find_cra(helix.end).chain
        start_auth_label = str(helix.start.res_id.seqid.num) + helix.start.res_id.seqid.icode
        end_auth_label = str(helix.end.res_id.seqid.num) + helix.end.res_id.seqid.icode
        start_id = chain[start_auth_label][0].label_seq
        end_id = end_chain[end_auth_label][0].label_seq
        left = min(start_id, end_id)
        right = max(start_id, end_id)
        if chain != end_chain:
            print('Helix ' + str(index) + ' in protein ' + id\
                               + ' is ill-defined. Unable to extract random coils.')
            return data
        else:
            secondary_structures.append((chain.name, left, right))
    
    for sheet in struct.sheets:
        for strand in sheet.strands:
            chain = struct[0].find_cra(strand.start).chain
            end_chain = struct[0].find_cra(strand.end).chain
            start_auth_label = str(strand.start.res_id.seqid.num) + strand.start.res_id.seqid.icode
            end_auth_label = str(strand.end.res_id.seqid.num) + strand.end.res_id.seqid.icode
            start_id = chain[start_auth_label][0].label_seq
            end_id = end_chain[end_auth_label][0].label_seq
            left = min(start_id, end_id)
            right = max(start_id, end_id)
            if chain != end_chain:
                print('Strand ' + strand.name + ' in sheet ' + sheet.name + ' in protein ' + id\
                                + ' is ill-defined. Unable to extract random coils.')
                return data
            else:
                secondary_structures.append((chain.name, left, right))
    
    # We want to scan through the secondary structures in order of chain, then by the starting sequence id.
    secondary_structures.sort(key=lambda x : (len(x[0]), x[0], x[1], x[2]))
    coil_id = 1
    ss_index = 0 # secondary structure index
    for chain in sequence.chain_start_indices:
        chain_object = struct[0].find_chain(chain)
        # If chain_object is None (which may happen if the whole chain is experimentally unconfirmed)
        if not chain_object:
            span = []
            chain_string = ""
        else:
            span = list(chain_object.get_polymer().first_conformer())
            chain_string = chain_object.get_polymer().make_one_letter_sequence()

        # Keep scanning through chain, iterating through helices and sheets/strands in the chain
        # until we run out of helices and sheets/strands that are in the chain.
        # We add any random coils that are before the current helix or strand pointed to by ss_index,
        # or add the rest of the random coils if there are no more helices or strands.
        while(True):
            # if ss_index points to the first helix or strand in the chain
            if ss_index == 0 or secondary_structures[ss_index - 1][0] != chain:
                coil_start = sequence.get_chain_start_id(chain)
            else:
                # coil starts after any other helix or strand that has come before the current helix or strand.
                coil_start = max(coil_start, secondary_structures[ss_index - 1][2] + 1)
            # if ss_index points to the first helix or strand in a subsequent chain,
            # or if we are out of helices or strands
            if ss_index >= len(secondary_structures) or secondary_structures[ss_index][0] != chain:
                coil_end = sequence.get_chain_end_id(chain)
            else:
                coil_end = secondary_structures[ss_index][1] - 1
            
            # In the case where the current helix/strand is at the start of a chain,
            # or there is a previous helix/strand at the end of the chain,
            # we get that coil_start = coil_end + 1, so we use this condition to ignore those cases.
            if coil_start <= coil_end:
                coil_sequence = sequence.get_chain_subsequence(chain, coil_start, coil_end)[0]
                annotated_sequence = sequence.get_chain_annotated_subsequence(span, chain_string, coil_start, coil_end)
                length = len(coil_sequence)
                unconfirmed = sequence.contains_unconfirmed_residues(chain, coil_start, coil_end)
                data.append((id, coil_id, chain, unconfirmed, coil_sequence, annotated_sequence,
                             coil_start, coil_end, length))
                coil_id += 1

            # If the current helix or strand is in the chain, then there may be more coils after
            # the current helix or strand, so we iterate to the next secondary structure.
            if ss_index < len(secondary_structures) and secondary_structures[ss_index][0] == chain:
                ss_index += 1
            else:
                break
    
    return data
