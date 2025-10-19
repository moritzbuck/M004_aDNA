import os, sys, shutil
from subprocess import call
import re
from tqdm import tqdm
import gzip
from ete3 import NCBITaxa
from Bio import Entrez
from speedict import Rdict, Options, AccessType
from random import choice

accession2taxid_path = "/home/moritz/scratch/dbs/nt/"

if os.path.exists(f"{accession2taxid_path}/all_accession2taxid.sdb"):
    db = Rdict(f"{accession2taxid_path}/all_accession2taxid.sdb", access_type=AccessType.read_only())

    
else:
    db = Rdict(f"{accession2taxid_path}/all_accession2taxid.sdb")
    with gzip.open(f"{accession2taxid_path}/accession2taxid/all_accession2taxid", "rt") as handle:
        handle.readline()
        for l in tqdm(handle, total=2_739_361_121):
            if not l.startswith("accession"):
                id_ = l.split()[0]
                taxid = l.split()[2]
                db[id_] = taxid
#speeDict

db_path = "/home/moritz/scratch/dbs/nt/bowtie2/"

db_files = os.listdir(db_path)
db_indexes = {l.split(".")[1] for l in db_files}
 

ncbi = NCBITaxa()
db_path = "/home/moritz/scratch/dbs/nt/bowtie2/"

db_files = os.listdir(db_path)
db_indexes = {l.split(".")[1] for l in db_files}


seq2tax = {}

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
    
def process_sam_file(sam):
    with open(sam) as handle:
        parsed_sam = [parse_samline(l) for l in tqdm(handle) if not l.startswith("@")]
    return parsed_sam

def match_counts(md):
    str_ = md.split(":")[-1]
    assert all([s in "^ATGC012345678N9" for s in str_])
    matches = sum([int(s) for s in re.split("[ATCGN^]", str_) if s])
    return matches

def parse_samline(l):
    saml = l.split()
    target = saml[2]
    read = saml[0]
    read_l = len(saml[9])
    md = [ll for ll in saml if ll.startswith("MD:")][0]
    matches = match_counts(md)
    return {'read'    : read,
            'target'  : target,
            'identity': matches/read_l}

def processing_sam(parsed_sam, identity) :
    parsed_sam = [k for k in parsed_sam if k['identity'] > identity]
    for p in tqdm(parsed_sam):
        taxid = db.get(str(p['target'].split(".")[0]), None)
        if taxid :
            if taxid not in seq2tax:
                seq2tax[taxid] = taxid2str(taxid)
            taxstrs = seq2tax[taxid]
            p['taxid'] = taxid
            p['taxstr_sci'] = taxstrs[0]
            p['taxstr_com'] = taxstrs[1]
    return parsed_sam
    
def run_subholi(fastq, idx, parse=True, identity = 0.92, preset = "new_holi"):
    tmp_fold = f"/scratch/moritz/{os.path.basename(fastq)}"

    if preset == "very_slow":
        # from the bear paper
        setting_string = "-t -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -k 10 "
    elif preset == "less_slow":
        # without the mismatch
        setting_string = "-t -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -k 10"
    elif preset == "new_holi":
        setting_string = "-t -k 1000"
    else :
        # bowtie2 default
        setting_string = ""

    bowtie_line  = f"bowtie2 -x {db_path}/nt.{idx} -U {fastq} {setting_string} --end-to-end --no-unal  --threads 20  -S {tmp_fold}/{idx}.sam"
    if os.path.exists(f"{tmp_fold}/{idx}.sam") :
        print(f"{idx} alreads done")
    else :
        print(f"doing {idx}")
        call(bowtie_line, shell = True)

    if parse:
        print(f"parsing {idx}")
        parsed_sam = process_sam_file(f"{tmp_fold}/{idx}.sam")
        return processing_sam(parsed_sam, identity)
    else :
        return True 
    
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

def lca(hit):
    tax_strs = [ db[h.reference_name.split(".")[0]] for h in hit if h.reference_name.split(".")[0] in db ]
    tax_strs = [taxid2str(h)[0] for h in tax_strs]
    taxes = [t.split(";") for t in tax_strs if t]
    inc = 1/len(hit)
    t_freqs = {}
    for i,k in enumerate(kranks):
        t_freqs[k] = {t[i] : 0 for t in taxes}
        for t in taxes:
            t_freqs[k][t[i]] += inc
    return {k : "" if(len(v) == 0) else max(v.items(), key = lambda x: x[1]) for k,v in t_freqs.items()}
            
            
def reholi(fastq):
    tmp_fold = f"/scratch/moritz/{os.path.basename(fastq)}"
    assert not os.path.exists(tmp_fold)
    os.makedirs(tmp_fold)
    parsed = []
    for idx in db_indexes:
        parsed += run_subholi(fastq,idx)

    hit_counts = {p['read'] : [] for p in tqdm(parsed)}
    for i,p in tqdm(enumerate(parsed)):
        if p['read'] not in hit_counts:
            hit_counts[p['read']] = []
        hit_counts[p['read']] += [i]

    lcas = {k :  lca(v)  for k,v  in tqdm(hit_counts.items())}
    consensus = [ ";".join([ v[k][0] if v[k][1] > 0.7 else ""  for k in sranks]) for v in tqdm(lcas.values())]
    
    with gzip.open(fastq,"rt") as handle:
        read_counts = len([l for l in handle if l.startswith("@")])

    tax_counts = {p : 0 for p in consensus}
    for p in consensus:
        tax_counts[p] +=1

def get_id(ali):
    md = ali.get_tag("MD")
    matches = match_counts(md)
    return matches/ali.query_length 

def post_proc(fastq):
    base = f"/home/moritz/scratch/{os.path.basename(fastq)}"
    sams = [f"{base}/{s}" for s in os.listdir(base) if s.endswith(".sam")]
    if len(sams) != 198 and not os.path.exists(f"{base}/final.sorted.bam"):
        print("not done mapping")
        return None
    elif os.path.exists(f"{base}/final.sorted.bam") :
        print("already mapped, already bam-ed")
        return 2
    else :
        print("done mapping, doing the bam-ing")
        pysam.merge("-o",f"{base}/final.bam", "-@", "20", *sams)
        pysam.sort("-o",f"{base}/final.sorted.bam", "-@", "20", f"{base}/final.bam")
        pysam.index("-@", "20", f"{base}/final.sorted.bam")
        print("done bam-ing")
        call(["rm", *sams])
        return 1
    
def sam2lca(fastq, id_cutoff = 0.8):
    if os.path.exists(f"/home/moritz/scratch/{os.path.basename(fastq)}/lca.best.json") :
        print("parsing already done, just loading from dumps")
        with open(f"/home/moritz/scratch/{os.path.basename(fastq)}/lca.best.json") as handle:
            tax_counts_best = json.load(handle)
        with open(f"/home/moritz/scratch/{os.path.basename(fastq)}/lca.consensus.json") as handle:
            tax_counts = json.load(handle)
        return  {'best_map' : tax_counts_best, 'consens_map' : tax_counts}
    with pysam.AlignmentFile(f"/home/moritz/scratch/{os.path.basename(fastq)}/final.sorted.bam") as handle:
        tot_len = sum([1 for l in handle])
    samfile = pysam.AlignmentFile(f"/home/moritz/scratch/{os.path.basename(fastq)}/final.sorted.bam")
    ali = next(samfile)
    current_read = ali.qname
    hits = [ali]
    data = {}
    for ali in tqdm(samfile, total = tot_len):
        if ali.qname != current_read:
            ids = [get_id(h) for h in hits]
            max_id = max(ids)
            bests = [h for h,id_ in zip(hits, ids) if id_ == max_id]
            if max_id >= id_cutoff and len(hits) >0:
                lca_best = lca(bests)
                lca_all = lca([h for id_,h in zip(ids, hits) if id_ >= id_cutoff])
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
    lca_best = lca(bests)
    lca_all = lca([h for id_,h in zip(ids, hits) if id_ > id_cutoff])
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

    with gzip.open(fastq,"rt") as handle:
        read_counts = len([l for l in handle if l.startswith("@")])

    tax_counts = {p : 0 for p in consensus}
    for p in consensus:
        tax_counts[p] +=1

    consensus_best = [ ";".join([ v['lca_best'][k][0] if v['lca_best'][k][1] > 0.9 else ""  for k in sranks[:-1]]) for v in tqdm(data.values())]

    tax_counts_best = {p : 0 for p in consensus_best}
    for p in consensus_best:
        tax_counts_best[p] +=1
    with open(f"/home/moritz/scratch/{os.path.basename(fastq)}/lca.best.json", "w") as handle:
        json.dump(tax_counts_best, handle, indent = 2)
    with open(f"/home/moritz/scratch/{os.path.basename(fastq)}/lca.consensus.json", "w") as handle:
        json.dump(tax_counts, handle, indent = 2)
    with open(f"/home/moritz/scratch/{os.path.basename(fastq)}/read_annotation.json", "w") as handle:
        json.dump(data, handle, indent = 2)
    samfile.close()
    return {'best_map' : tax_counts_best, 'consens_map' : tax_counts}

        

            
def merge_and_clean():
    tmp_fold = f"/scratch/moritz/{os.path.basename(fastq)}"
    for sam in os.listdir(tmp_fold):
        if sam.endswith(".sam") and not os.path.exists(f"{tmp_fold}/{sam.replace('.sam', '.bam')}"):
            call(f"""
        samtools view -bS {tmp_fold}/{sam} > {tmp_fold}/{sam.replace('.sam','.bam')} 
        """, shell=True)

    for bam in os.listdir(tmp_fold):
        if bam.endswith(".bam") and not os.path.exists(f"{tmp_fold}/{bam.replace('.bam', '.sbam')}"):
            call(f"""
            samtools sort -n -@ 22 -m 4G -o {tmp_fold}/{bam.replace('.bam', '.sbam')} {tmp_fold}/{bam}
        """, shell=True)

    call(f"samtools merge -@ 20 {tmp_fold}/*.sbam -o {tmp_fold}/final.bam", shell=True)
    call(f"filterBAM reassign --bam {tmp_fold}/final.bam -t 20 -i 0 -A 92 -M 30G -m 5G -n 10 -s 0.75 -o {tmp_fold}/final.reassign.bam", shell=True)
            

def do_all_m004_maps():
    raw_root = "/home/moritz/scratch/M004_aDNA/raws"
    fastqs = [f"{raw_root}/{f}" for f in  os.listdir(raw_root) if f.endswith(".highcomplex.fastq.gz") and not os.path.exists(f"/scratch/moritz/{os.path.basename(f)}/final.sorted.bam" )]
    for fastq in fastqs:
        tmp_fold = f"/scratch/moritz/{os.path.basename(fastq)}"
        os.makedirs(tmp_fold, exist_ok = True)
        parsed = []
        for idx in db_indexes:
            parsed += [run_subholi(fastq,idx, parse = False)]
 
def do_all_m004_post():
    raw_root = "/home/moritz/scratch/M004_aDNA/raws"
    fastqs = [f"{raw_root}/{f}" for f in  os.listdir(raw_root) if f.endswith(".highcomplex.fastq.gz")]
    data = {}
    for fastq in fastqs:
        print(f"checking if {fastq} is done mapping")
        if post_proc(fastq):
            print(f"parsing {fastq}")
            samp = os.path.basename(fastq).split(".highc")[0]
            data[samp] = sam2lca(fastq)



            
def classic_holi():
    tmp_fold = f"/scratch/moritz/{os.path.basename(fastq)}"
    assert not os.path.exists(tmp_fold)
    os.makedirs(tmp_fold)
    parsed = []
    for idx in db_indexes:
        parsed += run_subholi(fastq,idx, parse = False)
    assert all(parsed)
    merge_and_clean(fastq)




