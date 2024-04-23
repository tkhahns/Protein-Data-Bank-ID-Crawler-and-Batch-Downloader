import sqlite3
import os
import re
import extract
import attr
database = "./records/pdb_database_records.db"
rootdir = "./database"

def init(cur):
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.main_table + "(entry_id VARCHAR(5),\
                complex_type VARCHAR(25), contains_nonpolymer BOOL,\
                contains_antibody BOOL, chains VARCHAR, space_group VARCHAR(20),\
                Z_value INT, a FLOAT, b FLOAT, c FLOAT,\
                alpha FLOAT, beta FLOAT, gamma FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.entity_table + "(entry_id VARCHAR(5),\
                entity_id VARCHAR(5), entity_name VARCHAR(200), entity_type VARCHAR(25),\
                polymer_type VARCHAR(25), subchains VARCHAR)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.subchain_table +
                "(entry_id VARCHAR(5), entity_id VARCHAR(5), subchain_name VARCHAR(5),\
                parent_chain_name VARCHAR(5), chain_sequence VARCHAR, length INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.helix_table +
                "(entry_id VARCHAR(5), chain VARCHAR(5), helix_sequence VARCHAR, length INT)")


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
