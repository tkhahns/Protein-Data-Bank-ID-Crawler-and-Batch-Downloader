from enum import Enum

# All the table schemas that get produced in the database.
# First component is table name, second is all the attributes.
main_table = ("main", "(entry_id VARCHAR(5), complex_type VARCHAR(25), source_organism VARCHAR(200), chains VARCHAR,\
              space_group VARCHAR(20), Z_value INT, a FLOAT, b FLOAT, c FLOAT, alpha FLOAT, beta FLOAT, gamma FLOAT)")

entity_table = ("entities", "(entry_id VARCHAR(5), entity_id VARCHAR(5), entity_name VARCHAR(200),\
                entity_type VARCHAR(25), polymer_type VARCHAR(25), subchains VARCHAR)")

subchain_table = ("subchains", "(entry_id VARCHAR(5), entity_id VARCHAR(5), subchain_id VARCHAR(5),\
                  chain_id VARCHAR(5), chain_sequence VARCHAR, length INT)")

chain_table = ("chains", "(entry_id VARCHAR(5), chain_id VARCHAR(5), subchains VARCHAR,\
               chain_sequence VARCHAR, length INT)")

helix_table = ("helices", "(entry_id VARCHAR(5), chain_id VARCHAR(5), helix_sequence VARCHAR, start_position INT,\
               end_position INT, length INT)")

sheet_table = ("sheets", "(entry_id VARCHAR(5), sheet_id VARCHAR(5), number_strands INT, sense_sequence VARCHAR)")

strand_table = ("strands", "(entry_id VARCHAR(5), sheet_id VARCHAR(5),\
                strand_id VARCHAR(5), chain_id VARCHAR(5), strand_sequence VARCHAR, length INT)")

table_schemas = [main_table, entity_table, subchain_table, chain_table, helix_table, sheet_table, strand_table]

# Possible types of a complex, based on their entities
ComplexType = Enum('ComplexType', ['Other', 'SingleProtein',
                 'NucleicAcid', 'ProteinNA', 'Saccharide',
                 'ProteinSaccharide', 'SaccharideNA', 'ProteinSaccharideNA',
                 'Proteinmer', 'ComplexProtein'], start=0)