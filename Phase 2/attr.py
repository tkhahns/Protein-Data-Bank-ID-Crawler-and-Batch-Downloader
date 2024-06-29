from enum import Enum

# All the table schemas that get produced in the database.
# First component is table name, second is all the attributes.
main_table = ("main", "(entry_id VARCHAR(5) NOT NULL, complex_type VARCHAR(25) NOT NULL,\
              source_organism VARCHAR(200), chains VARCHAR, space_group VARCHAR(20),\
              Z_value INT, a FLOAT, b FLOAT, c FLOAT, alpha FLOAT, beta FLOAT, gamma FLOAT,\
              PRIMARY KEY (entry_id))")

entity_table = ("entities", "(entry_id VARCHAR(5) NOT NULL, entity_id VARCHAR(5) NOT NULL,\
                entity_name VARCHAR(200), entity_type VARCHAR(25), polymer_type VARCHAR(25), subchains VARCHAR,\
                PRIMARY KEY (entry_id, entity_id), FOREIGN KEY (entry_id) REFERENCES main (entry_id))")

subchain_table = ("subchains", "(entry_id VARCHAR(5) NOT NULL, entity_id VARCHAR(5) NOT NULL,\
                  subchain_id VARCHAR(5) NOT NULL, chain_id VARCHAR(5) NOT NULL,\
                  subchain_sequence VARCHAR, start_position INT, end_position INT, length INT,\
                  PRIMARY KEY (entry_id, subchain_id), FOREIGN KEY (entry_id) REFERENCES main (entry_id),\
                  FOREIGN KEY (entity_id) REFERENCES entities (entity_id),\
                  FOREIGN KEY (chain_id) REFERENCES chains (chain_id))")

chain_table = ("chains", "(entry_id VARCHAR(5) NOT NULL, chain_id VARCHAR(5) NOT NULL, subchains VARCHAR,\
               chain_sequence VARCHAR, start_position INT, end_position INT, length INT,\
               PRIMARY KEY (entry_id, chain_id), FOREIGN KEY (entry_id) REFERENCES main (entry_id))")

helix_table = ("helices", "(entry_id VARCHAR(5) NOT NULL, chain_id VARCHAR(5) NOT NULL, helix_sequence VARCHAR, start_position INT,\
               end_position INT, length INT, PRIMARY KEY (entry_id, chain_id, start_position, end_position),\
               FOREIGN KEY (entry_id) REFERENCES main (entry_id),\
               FOREIGN KEY (chain_id) REFERENCES chains (chain_id))")

sheet_table = ("sheets", "(entry_id VARCHAR(5) NOT NULL, sheet_id VARCHAR(5) NOT NULL, number_strands INT, sense_sequence VARCHAR,\
               PRIMARY KEY (entry_id, sheet_id), FOREIGN KEY (entry_id) REFERENCES main (entry_id))")

strand_table = ("strands", "(entry_id VARCHAR(5) NOT NULL, sheet_id VARCHAR(5) NOT NULL,\
                strand_id VARCHAR(5) NOT NULL, chain_id VARCHAR(5) NOT NULL,\
                strand_sequence VARCHAR, start_position INT, end_position INT, length INT,\
                PRIMARY KEY (entry_id, sheet_id, strand_id), FOREIGN KEY (entry_id) REFERENCES main (entry_id),\
                FOREIGN KEY (sheet_id) REFERENCES sheets (sheet_id)), FOREIGN KEY (chain_id) REFERENCES chains (chain_id)")

table_schemas = [main_table, entity_table, subchain_table, chain_table, helix_table, sheet_table, strand_table]

# Possible types of a complex, based on their entities
ComplexType = Enum('ComplexType', ['Other', 'SingleProtein',
                 'NucleicAcid', 'ProteinNA', 'Saccharide',
                 'ProteinSaccharide', 'SaccharideNA', 'ProteinSaccharideNA',
                 'Proteinmer', 'ComplexProtein'], start=0)