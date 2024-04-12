 # Protein Data Bank Analysis Code

 Project commissioned by Dr Wenqian Chen from the Department of Pharmacy and Pharmaceutical Sciences at the National University of Singapore as part of a lab assistant internship. The purpose of this project is to write tools to assist with managing and analysing all .pdb protein files from the RCSB Protein Data Bank (PDB). The project runs in 3 phases. Phase 1 is to develop a protocol for downloading and keeping and up-to-date copy of all the .pdb protein files from PDB. Phase 2 is to extract all the relevant information from the .pdb files. Phase 3 is to analyse all the relevant information.
 
## Phase 1

Previously, I used the PDB api and a batch download script to retrieve all the protein names and download their files. This has now been superceded by using `rsync` as recommended by PDB themselves, as it provides a simple interface for downloading and maintaining a local copy of their database. Instructions for using `rsync` with the database can be found [here](https://www.wwpdb.org/ftp/pdb-ftp-sites) The old phase 1 files are still here for archival purposes.

 The software uses the RCSB PDB search API to scrape a list of all protein IDs off of the data bank and writes it to `list_file.txt`. One can then run the bash script (sourced from PDB) to download all the PDB files; run `./batch_download.sh -f list_file.txt -p` in bash in the root directory, and a `.pdb.gz` archive file containing a `.ent` file of each protein will be downloaded to the `database` directory.

 The code uses the `node-fetch` library and is written as a JavaScript ES module. If running from source code, make sure to install `node-fetch` by running `npm install node-fetch` in the root directory. The code was tested on a Windows-x64 machine with Node 20.

 ## Phase 2

 The current plan is to use Python and SQLite3 to extract the relevant information from the .pdb files (id, name, chain structures and crystal structures, etc.) and store them in various tables in an SQL database.