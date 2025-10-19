import pysam
import argparse
import json
import sys
from tqdm import tqdm
import re
from speedict import Rdict, AccessType
from ete3 import NCBITaxa
from random import choice

def get_id(ali):
    md = ali.get_tag("MD")
    matches = match_counts(md)
    return matches/ali.query_length 

def match_counts(md):
    str_ = md.split(":")[-1]
    assert all([s in "^ATGC012345678N9" for s in str_])
    matches = sum([int(s) for s in re.split("[ATCGN^]", str_) if s])
    return matches

def lca_pair(a,b):
    diff = [i for i,aabb in enumerate(zip(a,b)) if aabb[0] != aabb[1]]
    if len(diff):
        return a[:diff[0]]
    else :
        return a

    lca = taxs[0]
    for tax in taxs[1:]:
        lca = lca_pair(lca, tax)
    lca = [l if l ==l else "" for l in lca]
    clust['lca_tax'] = ";".join(lca)

def lca(hit, db):
    tax_strs = [ db[h.reference_name.split(".")[0]] for h in hit if h.reference_name.split(".")[0] in db ]
    tax_strs = [taxid2str(h)[0] for h in tax_strs]
    taxes = [t.split(";") for t in tax_strs if t]
    inc = 1/len(hit)
    t_freqs = {}
    for i,k in enumerate(kranks):
        t_freqs[k] = {t[i] : 0 for t in taxes}
        for t in taxes:
            t_freqs[k][t[i]] += inc
    return {k : ('', 0.0) if(len(v) == 0) else max(v.items(), key = lambda x: x[1]) for k,v in t_freqs.items()}


kranks = ['superkingdom',
'kingdom',
'phylum',
'subphylum',
'class',
'subclass',
'infraclass',
'cohort',
'order',
'suborder',
'infraorder',
'superfamily',
'family',
'subfamily',
'genus',
'species']

sranks = ['superkingdom',
'phylum',
'class',
'order',
'family',
'genus',
'species']

ncbi = NCBITaxa()

_taxid2str = {}
def taxid2str(taxid):
    taxid = int(taxid)
    if taxid not in _taxid2str:
        if ncbi.get_taxid_translator([taxid]):
            lineage = ncbi.get_lineage(taxid)
            taxa = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            common_names = ncbi.get_common_names(lineage)
            common_name = ";".join([common_names.get(l, "") for l in lineage if ranks[l] in kranks])
            ranks = {v : k for k,v in ranks.items() if v in kranks}
            strs = [ "" if k not in ranks else taxa[ranks[k]] for k in kranks]
            taxstr = ";".join(strs)
            _taxid2str[taxid] = (taxstr, common_name)
        else :
            _taxid2str[taxid] = (None, None)
    return _taxid2str[taxid]



def main(args) :
    bam = args.bam
    id_cutoff = args.min_identity
    accession2taxid = args.accession2taxid
    prefix = args.prefix

    print(f"loading taxdb {accession2taxid}", sys.stderr)
    db = Rdict(f"{accession2taxid}", access_type=AccessType.read_only())


    with pysam.AlignmentFile(f"{bam}") as handle:
        tot_len = sum([1 for l in handle])
    samfile = pysam.AlignmentFile(f"{bam}")  
    current_read = None
    data = {}
    for ali in tqdm(samfile, total = tot_len):
        if current_read is None or ali.qname != current_read: 
            if current_read is None:
                    hits = [ali]
                    data = {}
                    current_read = ali.qname
            ids = [get_id(h) for h in hits]
            max_id = max(ids)
            bests = [h for h,id_ in zip(hits, ids) if id_ == max_id]
            if max_id >= id_cutoff and len(hits) >0:
                lca_best = lca(bests, db)
                lca_all = lca([h for id_,h in zip(ids, hits) if id_ >= id_cutoff], db)
                keep_one = choice(bests)
                data[current_read] = {
                    'lca_best' : lca_best,
                    'lca_all' : lca_all,
                    'nb_hits' : sum([id_ > id_cutoff for id_ in ids]),
                    'read_id' : keep_one.qname,
                    'read_seq' : keep_one.query,
                    'read_qual' : keep_one.qqual,
                    'reference_id' : keep_one.reference_name,
                    'reference_start' : keep_one.reference_start,
                    'reference_end' : keep_one.reference_end 
                }

            current_read = ali.qname
            hits = [ali]
        else :
            hits += [ali]
            lca_best = lca(bests, db)
            passers = [h for id_,h in zip(ids, hits) if id_ > id_cutoff]
            if len(passers) > 0:
                lca_all = lca(passers, db)
            else :
                lca_all = {k : ("", 1) for k in kranks}
            keep_one = choice(bests)            
            data[current_read] = {
                'lca_best' : lca_best,
                'lca_all' : lca_all,
                'nb_hits' : sum([id_ > id_cutoff for id_ in ids]),
                'read_id' : keep_one.qname,
                'read_seq' : keep_one.query,
                'read_qual' : keep_one.qqual,
                'reference_id' : keep_one.reference_name,
                'reference_start' : keep_one.reference_start,
                'reference_end' : keep_one.reference_end 
            }
    consensus = [ ";".join([ v['lca_all'][k][0] if v['lca_all'][k][1] > 0.9 else ""  for k in sranks[:-1]]) for v in tqdm(data.values())]

    tax_counts = {p : 0 for p in consensus}
    for p in consensus:
        tax_counts[p] +=1

    consensus_best = [ ";".join([ v['lca_best'][k][0] if v['lca_best'][k][1] > 0.9 else ""  for k in sranks[:-1]]) for v in tqdm(data.values())]

    tax_counts_best = {p : 0 for p in consensus_best}
    for p in consensus_best:
        tax_counts_best[p] +=1
    with open(f"{prefix}.lca.best.json", "w") as handle:
        json.dump(tax_counts_best, handle, indent = 2)
    with open(f"{prefix}.lca.consensus.json", "w") as handle:
        json.dump(tax_counts, handle, indent = 2)
    with open(f"{prefix}.read_annotation.json", "w") as handle:
        json.dump(data, handle, indent = 2)
    samfile.close()


parser = argparse.ArgumentParser(description="Assign LCA from BAM file")
parser.add_argument("--bam", required=True, help="Input BAM file path")
parser.add_argument("--accession2taxid", required=True, help="Path to accession2taxid Rdb file")
parser.add_argument("--min-identity", type=float, required=True, help="Minimum identity cutoff for filtering alignments")
parser.add_argument("--prefix", required=True, help="Output file prefix")

args = parser.parse_args()

main(args)
