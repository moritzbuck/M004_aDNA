import argparse
import json
import sys
from tqdm import tqdm
import re
from random import choice
import os
from speedict import Rdict, AccessType

from Bio import Entrez

# Set your email (required by NCBI for API requests)
Entrez.email = "your_email@example.com"

def get_taxid(gi):
    dbid = Entrez.read(Entrez.esearch(db="assembly", term=gi, idtype="acc", retmax = 5))['IdList'][0]
    md = Entrez.read(Entrez.esummary(db = "assembly", id = dbid))
    return md['DocumentSummarySet']['DocumentSummary'][0]['SpeciesTaxid']


def main(args) :
    dir_ = args.dir[0]
    accession2taxid = args.accession2taxid[0]
    ouput = args.output[0]
    db_dir = args.db_dir[0]
    print(f"finding index dbs in {db_dir}", sys.stderr)
    existing_bt_files = os.listdir(db_dir)
    # read the Rdict database into a dictionary 
    print(f"loading taxdb {accession2taxid}", sys.stderr)
    db = Rdict(f"{accession2taxid}")
    print(f"finding genome dbs in {dir_}", sys.stderr)
    genomes = os.listdir(dir_)

    has_been_indexed = lambda x : any([f.startswith(x) for f in existing_bt_files])
    to_process = []
    for genome in genomes :
        if not has_been_indexed(genome):
            print(f"PArsing {genome} and adding to taxdb")
            taxid = get_taxid(genome)
            fna = [f for f in os.listdir(f"{dir_}/{genome}") if f.endswith("_genomic.fna")]
            assert len(fna) == 1
            fna = fna[0]
            with open(f"{dir_}/{genome}/{fna}") as handle:
                for l in tqdm(handle):
                    if l[0] == ">":
                        db[l.split()[0].rstrip(">")] = taxid
            to_process += [genome]
        else :
            print(f"{genome} has been indexed allready")
    with open(ouput, "w") as handle:
        handle.write("\n".join(to_process) + "\n")



parser = argparse.ArgumentParser(description="searches thourgh the genome folder and adds the entries to the tax-db and checks which need to be indexed")
parser.add_argument("--dir", nargs=1, help="folder(s) with genomes in them")
parser.add_argument("--accession2taxid", nargs=1, help="Rdict database withj accession 2 taxids")
parser.add_argument("--output", nargs=1, help="output file")
parser.add_argument("--db_dir", nargs=1, help="folder with bt2 index files")
args = parser.parse_args()

main(args)
