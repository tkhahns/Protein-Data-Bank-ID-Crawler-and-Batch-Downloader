import sqlite3
import os
import re
import commands
from tqdm import tqdm
sql_database = "./Phase 2/records/pdb_database_records.db" # Location of output SQL database
rootdir = "./Phase 2/database" # Root directory of all the pdb files
verbose = False

if __name__ == "__main__":
    con = sqlite3.connect(sql_database)
    cur = con.cursor()
    commands.init_database(cur)

    for subdir, dirs, files in tqdm(os.walk(rootdir)):
        for file in files:
            path = os.path.join(subdir, file)
            if re.search('./*.cif.*', path):
                commands.check_file(cur, path, verbose=verbose)
        con.commit()

    con.close()
