import gemmi
from gemmi import cif
from typing import NamedTuple

class Monomer(NamedTuple):
    chain: str
    entity: int
    seq_id: int # An ordered labelling of the monomer/residue, not the same as its index in the polymer
    name: str
    expt_name: str # Same as name, but is labelled ? if experimentally unconfirmed
    hetero: str

class PolymerSequence:
    def __init__(self, doc: cif.Document):
        # We first extract all the relevant sequence-related info from the .cif file
        block = doc.sole_block()
        list_chain = list(block.find_loop("_pdbx_poly_seq_scheme.pdb_strand_id"))
        list_entity = list(block.find_loop("_pdbx_poly_seq_scheme.entity_id"))
        list_entity = [int(entity) for entity in list_entity]
        list_id = list(block.find_loop("_pdbx_poly_seq_scheme.seq_id"))
        list_id = [int(id) for id in list_id]
        list_monomer = list(block.find_loop("_pdbx_poly_seq_scheme.mon_id"))
        list_pdb_monomer = list(block.find_loop("_pdbx_poly_seq_scheme.pdb_mon_id"))
        list_hetero = list(block.find_loop("_pdbx_poly_seq_scheme.hetero"))
        
        temp_sequence = list(zip(list_chain, list_entity, list_id, list_monomer, list_pdb_monomer, list_hetero))
        # Convert list of tuples to list of named tuples for readability
        self.sequence = [Monomer(*monomer) for monomer in temp_sequence]
        self.bad_indices = []
        self.chain_start_indices = {}
        self.chain_end_indices = {}

        last_chain = ''
        index = 0
        while index < len(self.sequence):
            monomer = self.sequence[index]

            if monomer.hetero == 'y':
                while index + 1 < len(self.sequence) and self.sequence[index].seq_id == self.sequence[index + 1].seq_id:
                    self.sequence.pop(index + 1)
                    
            if monomer.expt_name == '?':
                self.bad_indices.append(index)

            if last_chain != monomer.chain:
                if last_chain != '':
                    self.chain_end_indices[last_chain] = index - 1
                last_chain = monomer.chain
                self.chain_start_indices[last_chain] = index
            index += 1
            
        self.chain_end_indices[last_chain] = index - 1
        self.chain_start_indices = dict(sorted(self.chain_start_indices.items(), key=lambda x : (len(x), x)))
        self.chain_end_indices = dict(sorted(self.chain_end_indices.items(), key=lambda x : (len(x), x)))

        monomer_sequence = [monomer.name for monomer in self.sequence]
        self.one_letter_code = sequence_3to1(monomer_sequence)
    
    def binary_search(self, left_index: int, right_index: int, target_label: int) -> int:
        """
        Helper method for finding the index of a desired residue in a residue sequence by means of binary search.
        This binary search is slightly modified to be faster than regular binary search.
        This modification only works if every value in the list we search through is unique,
        which is true when we let reduce all heterogeneities to single residues in our residue sequences.

        Keyword arguments:
        left_index -- index of the left bound of the sublist we're searching through
        right_index -- index of the right bound of the sublist we're searching through
        target_label -- the 'sequence id' of the residue we want to search for.
        """
        if right_index < left_index:
            raise Exception("Couldn't find index")
        centre_index = (left_index + right_index) // 2
        centre_label = self.sequence[centre_index].seq_id
        if centre_label == target_label:
            return centre_index
        
        difference = target_label - centre_label
        # By assuming we're searching over a list of distinct integers, we get that the target's index
        # is between centre_index and centre_index + difference
        next_index = centre_index + difference
        # Ensure that next_index are between the left_index and right_index
        next_index = max(left_index, next_index)
        next_index = min(right_index, next_index)
        next_label = self.sequence[next_index].seq_id

        if next_label == target_label:
            return next_index
        if target_label < centre_label:
            return self.binary_search(left_index, centre_index - 1, target_label)
        # centre_label < target_label
        return self.binary_search(centre_index + 1, right_index, target_label)
    
    def get_chain_sequence(self, chain: str) -> str:
        """
        Returns the full (unannotated) one-letter sequence of a given chain.
        """
        if chain not in self.chain_start_indices or chain not in self.chain_end_indices:
            return ''
        start = self.chain_start_indices[chain]
        end = self.chain_end_indices[chain]
        return self.one_letter_code[start:end + 1]
    
    def get_chain_subsequence(self, chain: str, start_id: int, end_id: int) -> tuple[str, int]:
        """
        Returns the (unannotated) one-letter sequence of
        a (possibly reversed) sublist of a given chain, and its signed length.
        """
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
    
    def get_chain_annotated_subsequence(self, span: list, chain_string: str, start_id: int, end_id: int) -> str:
        """
        Returns the annotated one-letter sequence of a sublist of a given chain.
        The full annotated one-letter sequence is derived from the gemmi function
        ResidueSpan.make_one_letter_sequence().
        """
        if len(span) == 0:
            return ""

        start_index = binary_search(span, 0, len(span) - 1, start_id)
        end_index = binary_search(span, 0, len(span) - 1, end_id)

        span_index_to_string_index = [i for i in range(len(chain_string)) if chain_string[i] != '-']
        string_start_index = span_index_to_string_index[start_index]
        string_end_index = span_index_to_string_index[end_index]
        if end_index >= start_index:
            return chain_string[string_start_index:string_end_index + 1]
        if end_index == 0:
            return self.one_letter_code[string_start_index::-1]
        return self.one_letter_code[string_start_index:string_end_index-1:-1]
    
    def contains_unconfirmed_residues(self, chain: str, start_id: int, end_id: int) -> int:
        """
        Checks if a sublist of a chain (between start_id and end_id) contains an experimentally unconfirmed residue.
        Note that this function returns 0 or 1 instead of a bool, as sqlite does not natively support booleans.
        This may change if we move to a different SQL engine.
        """
        chain_start = self.chain_start_indices[chain]
        chain_end = self.chain_end_indices[chain]
        
        for index in self.bad_indices:
            if chain_start <= index <= chain_end:
                return 1
        return 0
    
    def get_chain_start_id(self, chain: str) -> int:
        """
        Returns the sequence id of the starting residue of a chain (not its index).
        """
        return self.sequence[self.chain_start_indices[chain]].seq_id
    
    def get_chain_end_id(self, chain: str) -> int:
        """
        Returns the sequence id of the ending residue of a chain (not its index).
        """
        return self.sequence[self.chain_end_indices[chain]].seq_id
    
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
            return ("MULTIPLE CHAINS ERROR", 0)
        return self.get_chain_subsequence(chain.name, start_pos, end_pos)
    
def binary_search(span: list[gemmi.Residue], left_index: int, right_index: int, target_label: int) -> int:
    """
    Helper method for finding the index of a desired residue in a residue span by means of binary search.
    This binary search is slightly modified to be faster than regular binary search.
    This modification only works if every value in the list we search through is unique,
    which is true when we let reduce all heterogeneities to single residues in our residue spans.

    Keyword arguments:
    left_index -- index of the left bound of the sublist we're searching through
    right_index -- index of the right bound of the sublist we're searching through
    target_label -- the 'sequence id' of the residue we want to search for.
    """
    if right_index < left_index:
        raise Exception("Couldn't find index")
    if target_label < span[left_index].label_seq:
        return left_index
    if target_label > span[right_index].label_seq:
        return right_index
    centre_index = (left_index + right_index) // 2
    centre_label = span[centre_index].label_seq
    if centre_label == target_label:
        return centre_index
    
    difference = target_label - centre_label
    # By assuming we're searching over a list of distinct integers, we get that the target's index
    # is between centre_index and centre_index + difference
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
                'TYR': 'Y', 'VAL': 'V', 'DA': 'A',
                'DT': 'T', 'DG': 'G', 'DC': 'C'}

def letter_code_3to1(polymer: str) -> str:
    if (polymer in three_to_one.keys()):
        return three_to_one[polymer]
    return 'X'

def sequence_3to1(sequence: list[str]) -> str:
    return ''.join([letter_code_3to1(polymer) for polymer in sequence])
