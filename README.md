 # Protein Data Bank ID Crawler & Batch Downloader

 Project commissioned by Dr Wenqian Chen from the Department of Pharmacy and Pharmaceutical Sciences at the National University of Singapore. The purpose of this small project is to develop a simple(r) way to batch download all the protein .pdb files off of the Protein Data Bank (PDB).
 
 So far, this software uses the RCSB PDB search API to scrape a list of all protein IDs off of the data bank and writes it to `list_file.txt`. One can then run the bash script (sourced from PDB) to download all the PDB files; run `./batch_download.sh -f list_file.txt -p` in bash in the root directory, and a `.pdb.gz` archive file containing a `.ent` file of each protein will be downloaded to the `database` directory.

 ## Requirements

 The code uses the `node-fetch` library and is written as a JavaScript ES module. If running from source code, make sure to install `node-fetch` by running `npm install node-fetch` in the root directory. The code was tested on a Windows-x64 machine with Node 20.