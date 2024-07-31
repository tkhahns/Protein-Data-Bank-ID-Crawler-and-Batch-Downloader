import sqlite3
import gemmi
from gemmi import cif
from database import table_schemas

def init_database(cur: sqlite3.Cursor):
    for table_schema in table_schemas:
        cur.execute(table_schema.create_table())

def insert_file(cur: sqlite3.Cursor, file_path: str, verbose = True):
    #try:
    if verbose:
        print(file_path)
    struct = gemmi.read_structure(file_path)
    doc = cif.read(file_path)
    res = cur.execute("SELECT entry_id FROM " + table_schemas[0].name\
                        + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
    if not res.fetchone(): # if there is no row in the main table with such entry ID
        for table_scheme in table_schemas:
            for data in table_scheme.extract_data(struct, doc):
                statement = table_scheme.insert_row(data)
                cur.execute(statement, data)
    #except Exception as error:
    #    print(struct.name)
    #    print(error)

# def update_database()