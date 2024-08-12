# Protein Data Bank Analysis Code

 Project commissioned by Dr Wenqian Chen from the Department of Pharmacy and Pharmaceutical Sciences at the National University of Singapore as part of a lab assistant internship. The purpose of this project is to write tools to assist with managing and analysing all protein mmCIF files from the RCSB Protein Data Bank (PDB). The project runs in 3 phases. Phase 1 is to develop a protocol for downloading and keeping and up-to-date copy of all the mmCIF protein files from PDB. Phase 2 is to extract all the relevant information from the mmCIF files. Phase 3 is to analyse all the relevant information.

## Phase 1

 Previously, I used the PDB api and a batch download script to retrieve all the protein names and download their files. This has now been superseded by using `rsync` as recommended by PDB themselves, as it provides a simple interface for downloading and maintaining a local copy of their database. Instructions for using `rsync` with the database can be found [here](https://www.wwpdb.org/ftp/pdb-ftp-sites) The old phase 1 files are still here for archival purposes.

 The software uses the RCSB PDB search API to scrape a list of all protein IDs off of the data bank and writes it to `list_file.txt`. One can then run the bash script (sourced from PDB) to download all the PDB files; run `./batch_download.sh -f list_file.txt -p` in bash in the root directory, and a `.pdb.gz` archive file containing a `.ent` file of each protein will be downloaded to the `database` directory.

 The code uses the `node-fetch` library and is written as a JavaScript ES module. If running from source code, make sure to install `node-fetch` by running `npm install node-fetch` in the root directory. The code was tested on a Windows-x64 machine with Node 20.

## Phase 2

 We use Python and SQLite3 to extract the relevant information from the .pdb files (id, name, cell structure, primary chain structure, secondary alpha helix and beta sheet structures, component entities, etc.) and store them in various tables in an SQL database. If you wish to run this code yourself, make sure to change the `database` and `rootdir` variables in `main.py` before running `main.py` through Python. The GEMMI Python library is used to extract molecule structure information.

 See GEMMI documentation [here](https://gemmi.readthedocs.io/en/latest/index.html).

 I recommend using [this guide](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction) and [this dictionary resource](https://mmcif.wwpdb.org) to understand the .cif file structure.

### SQL table information

 The SQL tables produced have the following schema. **Bold** indicates a primary key, *italics* indicate foreign key, and ***both*** indicate foreign primary key.

 | Table name | Attributes |
 | ---------- | ---------- |
 | main       | **entry_id**, complex_type, source_organism, chains, space_group, Z_value, a, b, c, alpha, beta, gamma |
 | entities   | ***entry_id***, **entity_id**, entity_name, entity_type, polymer_type, subchains |
 | subchains  | ***entry_id***, *entity_id*, **subchain_id**, *chain_id*, chain_sequence, start_position, end_position, length |
 | chains     | ***entry_id***, **chain_id**, subchains, chain_sequence, start_position, end_position, length |
 | helices    | ***entry_id***, ***chain_id***, helix_sequence, **start_position**, **end_position**, length |
 | sheets     | ***entry_id***, **sheet_id**, number_strands, sense_sequence |
 | strands    | ***entry_id***, ***sheet_id***, **strand_id**, *chain_id*, strand_sequence, start_position, end_position, length |

 An explanation on how sequences work is warranted, despite how simple they may seem. All sequences (chain, subchain, helix or strand) consist of the one letter code of each amino acid residue of the chain/subchain/helix/strand span. Details about what each letter represents can be found [here](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.one_letter_code.html). Besides the Latin alphabet letters, there can also be dashes in the sequence, representing either a segment of the sequence that hasn't been experimentally confirmed, or a link between two independent components of the span.

 Polymers may contain microhomogeneities, i.e. alternative residues can occupy the same sequence ID without any major change in the properties of the polymer. In such cases, the sequence contains the 'first conformer'. The start and end positions of the sequence indicate the sequence ID that the start and end residues occupy in the span. Note that sequence IDs of a span doesn't count from 1 and up; it can start on any number and may skip numbers as the original author sees fit. The length counts how many residues are in the span, only counting the first conformer of a set of microhomogeneities and currently excluding experimentally unconfirmed residues. The length is negative for a helix or strand sequence if the helix or strand goes in the opposite direction of the parent chain.

 Each residue actually has two sequence IDs: the primary sequence ID and the author sequence ID. The primary sequence ID is used, as it runs in a strictly increasing order along the span, which makes it useful for coding with. There may be multiple residues in the span with the same numerical author sequence ID, and differ only in their icode (see the Python class property gemmi.SeqId.icode). Hence, the author sequence ID may not necessarily be increasing, which makes it harder to work with in code.

 The 'sense sequence' in the sheets table indicate whether adjacent strands in the sheet are running parallel or antiparallel with each other. So if a sheet has sense sequence "PPA", it means that the sheet contains four strands, with the first three strands being parallel with each other and the third and fourth strand being antiparallel with each other.
