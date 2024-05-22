import sqlite3
import os
import re
import extract
import attr
database = "./Phase 2/records/pdb_database_records.db" # Location of output SQL database
rootdir = "./Phase 2/database" # Root directory of all the pdb files

def init(cur):
    for table_schema in attr.table_schemas:
        cur.execute("CREATE TABLE IF NOT EXISTS " + table_schema[0] + table_schema[1])

if __name__ == "__main__":
    con = sqlite3.connect(database)
    cur = con.cursor()
    init(cur)

    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            path = os.path.join(subdir, file)
            if re.search('.*\.cif.*', path):
                extract.insert_into_all_tables(path, cur)

    con.commit()
    con.close()
