import sqlite3
from table import Table
from attributes import Attributes
import extract

# All the table schemas that get produced in the database.
# First component is table name, second is all the attributes.

entry_id = ("entry_id", "VARCHAR(5) NOT NULL")
chain_id = ("chain_id", "VARCHAR(5) NOT NULL")
sheet_id = ("sheet_id", "VARCHAR(5) NOT NULL")
unconfirmed = ("contains_experimentally_unconfirmed_residues", "INT")
start_id = ("start_id", "INT")
end_id = ("end_id", "INT")
length = ("length", "INT")

main_table_attributes = Attributes[extract.MainData]\
    ([entry_id, ("complex_type", "VARCHAR(25) NOT NULL"), ("structure_title", "VARCHAR"),
      ("source_organism", "VARCHAR(200)"), ("revision_date", "VARCHAR(50)"), ("chains", "VARCHAR"),
      ("space_group", "VARCHAR(20)"), ("Z_value", "INT"), ("a", "FLOAT"), ("b", "FLOAT"), ("c", "FLOAT"),
      ("alpha", "FLOAT"), ("beta", "FLOAT"), ("gamma", "FLOAT")],
      primary_keys=["entry_id"])
main_table = Table("main", main_table_attributes, extract.insert_into_main_table)

experimental_table_attributes = Attributes[extract.ExperimentalData]\
    ([entry_id, ("Matthews_coefficient", "FLOAT"), ("percent_solvent_content", "FLOAT"),
      ("crystal_growth_method", "VARCHAR"), ("crystal_growth_procedure", "VARCHAR"),
      ("crystal_growth_apparatus", "VARCHAR"), ("crystal_growth_atmosphere", "VARCHAR"),
      ("crystal_growth_pH", "FLOAT"), ("crystal_growth_temperature", "FLOAT")],
      primary_keys=["entry_id"],
      foreign_keys={"entry_id": ("main", "entry_id")})
experimental_table = Table("experimental", experimental_table_attributes, extract.insert_into_experimental_table)

entity_table_attributes = Attributes[extract.EntityData]\
    ([entry_id, ("entity_id", "VARCHAR(5) NOT NULL"), ("entity_name", "VARCHAR(200)"),
      ("entity_type", "VARCHAR(25)"), ("polymer_type", "VARCHAR(25)"), ("subchains", "VARCHAR")],
      primary_keys=["entry_id", "entity_id"],
      foreign_keys={"entry_id": ("main", "entry_id")})
entity_table = Table("entities", entity_table_attributes, extract.insert_into_entity_table)

chain_table_attributes = Attributes[extract.ChainData]\
    ([entry_id, chain_id, ("subchains", "VARCHAR"), unconfirmed, ("chain_sequence", "VARCHAR"),
      ("annotated_chain_sequence", "VARCHAR"), start_id, end_id, length, ("author_start_id", "INT"), ("author_end_id", "INT")],
      primary_keys=["entry_id", "chain_id"],
      foreign_keys={"entry_id": ("main", "entry_id")})
chain_table = Table("chains", chain_table_attributes, extract.insert_into_chain_table)

subchain_table_attributes = Attributes[extract.SubchainData]\
    ([entry_id, ("entity_id", "VARCHAR(5) NOT NULL"), ("subchain_id", "VARCHAR(5) NOT NULL"), chain_id,
      ("subchain_sequence", "VARCHAR"), ("annotated_subchain_sequence", "VARCHAR"),
      start_id, end_id, length],
      primary_keys=["entry_id", "subchain_id"],
      foreign_keys={"entry_id": ("main", "entry_id"), "entity_id": ("entities", "entity_id"),
                    "chain_id": ("chains", "chain_id")})
subchain_table = Table("subchains", subchain_table_attributes, extract.insert_into_subchain_table)

helix_table_attributes = Attributes[extract.HelixData]\
    ([entry_id, ("helix_id", "INT"), chain_id, ("helix_sequence", "VARCHAR"), ("helix_type", "VARCHAR"),
      start_id, end_id, length],
      primary_keys=["entry_id", "helix_id"],
      foreign_keys={"entry_id": ("main", "entry_id"), "chain_id": ("chains", "chain_id")})
helix_table = Table("helices", helix_table_attributes, extract.insert_into_helix_table)

sheet_table_attributes = Attributes[extract.SheetData]\
    ([entry_id, sheet_id, ("number_strands", "INT"), ("sense_sequence", "VARCHAR")],
     primary_keys=["entry_id", "sheet_id"],
     foreign_keys={"entry_id": ("main", "entry_id")})
sheet_table = Table("sheets", sheet_table_attributes, extract.insert_into_sheet_table)

strand_table_attributes = Attributes[extract.StrandData]\
    ([entry_id, sheet_id, ("strand_id", "VARCHAR(5) NOT NULL"), chain_id,
      ("strand_sequence", "VARCHAR"), start_id, end_id, length],
      primary_keys=["entry_id", "sheet_id", "strand_id"],
      foreign_keys={"entry_id": ("main", "entry_id"), "sheet_id": ("sheets", "sheet_id"),
                    "chain_id": ("chains", "chain_id")})
strand_table = Table("strands", strand_table_attributes, extract.insert_into_strand_table)

coil_table_attributes = Attributes[extract.CoilData]\
    ([entry_id, ("coil_id", "INT"), chain_id, unconfirmed, ("coil_sequence", "VARCHAR"),
      ("annotated_coil_sequence", "VARCHAR"), start_id, end_id, length],
      primary_keys=["entry_id", "coil_id"],
      foreign_keys={"entry_id": ("main", "entry_id"), "chain_id": ("chains", "chain_id")})
coil_table = Table("coils", coil_table_attributes, extract.insert_into_coil_table)

# Define secondary_structures table attributes
secondary_structures_table_attributes = Attributes[extract.HelixData]\
    ([entry_id,
      ("helix_id", "INT"),
      chain_id,
      ("helix_sequence", "VARCHAR"),
      ("helix_type", "VARCHAR"),
      start_id,
      end_id,
      length],
     primary_keys=["entry_id", "helix_id"],
     foreign_keys={
         "entry_id": ("main", "entry_id"),
         "chain_id": ("chains", "chain_id"),
         "helix_id": ("helices", "helix_id")
     })

# Create the secondary_structures table using the new attributes and an insert function in extract module.
secondary_structures_table = Table("secondary_structures",
                                   secondary_structures_table_attributes,
                                   extract.insert_into_secondary_structures_table)

table_schemas: list[Table] = [main_table, experimental_table, entity_table, chain_table,
                              subchain_table, helix_table, sheet_table, strand_table, coil_table,
                              secondary_structures_table]

def insert_into_table(cur: sqlite3.Cursor, table_name: str, data):
    """
    Inserts the given data into the given table
    """
    args = ', '.join(['?' for i in range(len(data))])
    query = f'INSERT INTO {table_name} VALUES({args})'
    cur.execute(query, data)

def retrieve_from_table(cur: sqlite3.Cursor, table_name: str, entry_id: str):
    """
    Retrieves all rows from a given table with the specified entry id.
    """
    result = cur.execute(f'SELECT * FROM {table_name} WHERE entry_id = {entry_id}')
    return result.fetchall()