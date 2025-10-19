import pysam
import argparse
import json
import sys
from tqdm import tqdm
import re
from random import choice
import pandas



def main(args) :
    bests = args.best if args.best else []
    consensuses = args.consensus if args.consensus else []
    annotations = args.annotation if args.annotation else []

    best_dict = {}
    consensus_dict = {}

    for f in bests :
        sample = f.split("/")[-1].split(".lca.best.")[0]
        with open(f) as handle:
            best_dict[sample] = json.load(handle)
    best_df = pandas.DataFrame.from_dict(best_dict).fillna(0)
    best_df.to_csv("lca.best.tsv", sep = "\t", index_label="Taxon")
    for f in consensuses :
        sample = f.split("/")[-1].split(".lca.consensus.")[0]
        with open(f) as handle:
            consensus_dict[sample] = json.load(handle)
    consensus_df = pandas.DataFrame.from_dict(consensus_dict).fillna(0)
    consensus_df.to_csv("lca.consensus.tsv", sep = "\t", index_label="Taxon")
    reads = []
    for f in annotations :
        sample = f.split("/")[-1].split(".read_annotation.")[0]
        with open(f) as handle:
            data = json.load(handle)
        for k,v in data.items():
            tax_str = ":".join([ v['lca_best'][k][0] if v['lca_best'][k][1] > 0.9 else ""  for k in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'][:-1]])
            read = f"@{k};{sample};{tax_str};{v['reference_id']}:{v['reference_start']}:{v['reference_end']}\n{v['read_seq']}\n+\n{v['read_qual']}\n"
            reads.append(read)
    with open("reads.fastq", "w") as handle:
        handle.writelines(reads)

parser = argparse.ArgumentParser(description="Assign LCA from BAM file")
parser.add_argument("--best", nargs="+", help="jsons with best LCA assignments")
parser.add_argument("--consensus", nargs="+", help="jsons with consensus LCA assignments")
parser.add_argument("--annotation", nargs="+", help="jsons with read annotations")

args = parser.parse_args()

main(args)
