import gemmi
from gemmi import cif

class PolymerSequence:
    def __init__(self, doc: cif.Document):
        block = doc.sole_block()
        list_chain = list(block.find_loop("_pdbx_poly_seq_scheme.pdb_strand_id"))
        list_entity = list(block.find_loop("_pdbx_poly_seq_scheme.entity_id"))
        list_entity = [int(entity) for entity in list_entity]
        list_id = list(block.find_loop("_pdbx_poly_seq_scheme.seq_id"))
        list_id = [int(id) for id in list_id]
        list_monomer = list(block.find_loop("_pdbx_poly_seq_scheme.mon_id"))
        list_hetero = list(block.find_loop("_pdbx_poly_seq_scheme.hetero"))
        
        self.sequence = list(zip(list_chain, list_entity, list_id, list_monomer, list_hetero))
        self.chain_start_indices = {}
        self.chain_end_indices = {}
        last_chain = ''
        index = 0
        while index < len(self.sequence):
            monomer = self.sequence[index]
            if monomer[4] == 'y':
                while self.sequence[index][2] == self.sequence[index + 1][2]:
                    self.sequence.pop(index + 1)
            if last_chain != monomer[0]:
                if last_chain != '':
                    self.chain_end_indices[last_chain] = index - 1
                last_chain = monomer[0]
                self.chain_start_indices[last_chain] = index
            index += 1
        self.chain_end_indices[last_chain] = index - 1

        amino_acid_sequence = [monomer[3] for monomer in self.sequence]
        self.one_letter_code = sequence_3to1(amino_acid_sequence)
    
    def binary_search(self, left_index: int, right_index: int, target_label: int) -> int:
        """
        Helper method for finding the index of a desired element in a residue span by means of binary search.
        """
        if right_index < left_index:
            raise Exception("Couldn't find index")
        centre_index = (left_index + right_index) // 2
        centre_label = self.sequence[centre_index][2]
        if centre_label == target_label:
            return centre_index
        difference = target_label - centre_label
        next_index = centre_index + difference
        # ensure that next_index are between the left_index and right_index
        next_index = max(left_index, next_index)
        next_index = min(right_index, next_index)
        next_label = self.sequence[next_index][2]
        if next_label == target_label:
            return next_index
        
        if target_label < centre_label:
            return self.binary_search(left_index, centre_index - 1, target_label)
        # centre_label < target_label
        return self.binary_search(centre_index + 1, right_index, target_label)
    
    def get_chain_sequence(self, chain: str) -> str:
        if chain not in self.chain_start_indices:
            return ''
        start = self.chain_start_indices[chain]
        end = self.chain_end_indices[chain]
        return self.one_letter_code[start:end + 1]
    
    def get_chain_subsequence(self, chain: str, start_id: int, end_id: int) -> tuple[str, int]:
        if chain not in self.chain_start_indices:
            return ''
        chain_start = self.chain_start_indices[chain]
        chain_end = self.chain_end_indices[chain]
        start_index = self.binary_search(chain_start, chain_end, start_id)
        end_index = self.binary_search(chain_start, chain_end, end_id)
        if end_index >= start_index:
            return self.one_letter_code[start_index:end_index + 1], end_index - start_index + 1
        if end_index == 0:
            return self.one_letter_code[start_index::-1], -start_index - 1
        return self.one_letter_code[start_index:end_index-1:-1], end_index - start_index - 1
    
    def get_helix_sequence(self, helix: gemmi.Helix, struct: gemmi.Structure) -> str:
        chain = struct[0].find_cra(helix.start).chain
        end_chain = struct[0].find_cra(helix.end).chain
        start_auth_label = str(helix.start.res_id.seqid.num) + helix.start.res_id.seqid.icode
        end_auth_label = str(helix.end.res_id.seqid.num) + helix.end.res_id.seqid.icode
        start_pos = chain[start_auth_label][0].label_seq
        end_pos = end_chain[end_auth_label][0].label_seq

        if (chain != end_chain):
            return "MULTIPLE CHAINS ERROR"
        return self.get_chain_subsequence(chain.name, start_pos, end_pos)[0]
    
    def get_strand_sequence(self, strand: gemmi.Sheet.Strand, struct: gemmi.Structure) -> tuple[str, int]:
        chain = struct[0].find_cra(strand.start).chain
        end_chain = struct[0].find_cra(strand.end).chain
        start_auth_label = str(strand.start.res_id.seqid.num) + strand.start.res_id.seqid.icode
        end_auth_label = str(strand.end.res_id.seqid.num) + strand.end.res_id.seqid.icode
        start_pos = chain[start_auth_label][0].label_seq
        end_pos = end_chain[end_auth_label][0].label_seq

        if (chain != end_chain):
            return ("MULTIPLE CHAINS ERROR", -1)
        return self.get_chain_subsequence(chain.name, start_pos, end_pos)
    
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
    return 'X'

def sequence_3to1(sequence: list[str], start: int = 0, end: int = -1) -> str:
    return ''.join([letter_code_3to1(polymer) for polymer in sequence[start:end]])
