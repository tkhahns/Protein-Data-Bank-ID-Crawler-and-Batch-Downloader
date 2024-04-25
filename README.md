 # Protein Data Bank Analysis Code

 Project commissioned by Dr Wenqian Chen from the Department of Pharmacy and Pharmaceutical Sciences at the National University of Singapore as part of a lab assistant internship. The purpose of this project is to write tools to assist with managing and analysing all protein mmCIF files from the RCSB Protein Data Bank (PDB). The project runs in 3 phases. Phase 1 is to develop a protocol for downloading and keeping and up-to-date copy of all the mmCIF protein files from PDB. Phase 2 is to extract all the relevant information from the mmCIF files. Phase 3 is to analyse all the relevant information.
 
## Phase 1

Previously, I used the PDB api and a batch download script to retrieve all the protein names and download their files. This has now been superceded by using `rsync` as recommended by PDB themselves, as it provides a simple interface for downloading and maintaining a local copy of their database. Instructions for using `rsync` with the database can be found [here](https://www.wwpdb.org/ftp/pdb-ftp-sites) The old phase 1 files are still here for archival purposes.

 The software uses the RCSB PDB search API to scrape a list of all protein IDs off of the data bank and writes it to `list_file.txt`. One can then run the bash script (sourced from PDB) to download all the PDB files; run `./batch_download.sh -f list_file.txt -p` in bash in the root directory, and a `.pdb.gz` archive file containing a `.ent` file of each protein will be downloaded to the `database` directory.

 The code uses the `node-fetch` library and is written as a JavaScript ES module. If running from source code, make sure to install `node-fetch` by running `npm install node-fetch` in the root directory. The code was tested on a Windows-x64 machine with Node 20.

 ## Phase 2

 We use Python and SQLite3 to extract the relevant information from the .pdb files (id, name, cell structure, primary chain structure, secondary alpha helix and beta sheet structures, component entities, etc.) and store them in various tables in an SQL database. If you wish to run this code yourself, make sure to change the `database` and `rootdir` variables in `main.py` before running `main.py` through Python. The GEMMI Python library is used to extract molecule structure information.

See GEMMI documentation [here](https://gemmi.readthedocs.io/en/latest/index.html).

I recommend using [this guide](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction) and [this dictionary resource](https://mmcif.wwpdb.org) to understand the .cif file structure.