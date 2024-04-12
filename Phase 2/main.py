import sqlite3
import os
import re
import extract
import attr
database = "./database/pdb_database.db"
rootdir = "./database"

def init(cur):
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.main_table + "(id VARCHAR(5),\
                protein_atoms_number INT, nucleic_acid_atoms_number INT,\
                space_group VARCHAR(20), Z_value INT, a FLOAT, b FLOAT, c FLOAT,\
                alpha FLOAT, beta FLOAT, gamma FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.entity_table + "(id VARCHAR(5),\
                entity_name VARCHAR(50), entity_type VARCHAR(50), chains VARCHAR, full_sequence VARCHAR)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.chain_table + "(id VARCHAR(5),\
                chain_name VARCHAR(5), chain_sequence VARCHAR)")
    cur.execute("CREATE TABLE IF NOT EXISTS " + attr.helix_table + "(id VARCHAR(5),\
                chain_name VARCHAR(5), helix_sequence VARCHAR)")


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
