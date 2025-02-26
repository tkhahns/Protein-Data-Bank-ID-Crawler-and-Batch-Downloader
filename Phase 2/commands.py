import sqlite3
import gemmi
from gemmi import cif
from database import table_schemas
from polymer_sequence import PolymerSequence

def init_database(cur: sqlite3.Cursor):
    for table_schema in table_schemas:
        cur.execute(table_schema.create_table())

def check_file(cur: sqlite3.Cursor, file_path: str, verbose: bool = True):
    try:
        if verbose:
            print("Checking " + file_path)
        struct = None
        struct = gemmi.read_structure(file_path)
        doc = cif.read(file_path)
        sequence = PolymerSequence(doc)

        # Check if protein file exists in database
        res = cur.execute("SELECT entry_id FROM " + table_schemas[0].name\
                            + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
        if not res.fetchone(): # if there is no row in the main table with such entry ID
            if verbose:
                print("Adding " + file_path)
            insert_file(cur, struct, doc, sequence)

        else: # Check if protein file data is up to date
            block = doc.sole_block()
            revision_date = block.find_value("_pdbx_audit_revision_history.revision_date")
            if revision_date is None:
                revision_date = block.find_loop("_pdbx_audit_revision_history.revision_date")[-1]
            res = cur.execute("SELECT revision_date FROM " + table_schemas[0].name\
                              + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
            if res.fetchone()[0] < revision_date:
                if verbose:
                    print("Updating " + file_path)
                update_file(cur, struct, doc, sequence)

            else: # Check that protein file data did not get corrupted
                res = cur.execute("SELECT entry_id FROM " + table_schemas[-1].name\
                                + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
                # if there is no row in the last table (coils) with such entry ID, then something went wrong.
                # I checked and every protein has some rows in the coils table.
                if not res.fetchone():
                    if verbose:
                        print("Data corrupted, fixing " + file_path)
                    update_file(cur, struct, doc, sequence)

    except Exception as error:
        if struct is not None:  
            print(struct.name)
            print(error)
        else:
            print(error)
            

def insert_file(cur: sqlite3.Cursor, struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence):
    for table_scheme in table_schemas:
        for data in table_scheme.extract_data(struct, doc, sequence):
            statement = table_scheme.insert_row(data)
            cur.execute(statement, data)

def update_file(cur: sqlite3.Cursor, struct: gemmi.Structure, doc: cif.Document, sequence: PolymerSequence):
    """
    Used to add data to all tables if the given protein only has data in some tables, or if data is not up to date.
    This may happen if regular file insertion was interrupted.
    """
    for table_scheme in table_schemas:
        cur.execute("DELETE FROM " + table_scheme.name + " WHERE entry_id = '" + struct.info["_entry.id"] + "'")
        for data in table_scheme.extract_data(struct, doc, sequence):
            statement = table_scheme.insert_row(data)
            cur.execute(statement, data)