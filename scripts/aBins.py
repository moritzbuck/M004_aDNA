from Bio import SeqIO
import os, sys, shutil
import pandas
from numpy import logical_and

pydam = pandas.read_csv("../Ancient_DNA/pydamage/filter/MEGAHIT-coassembly-bbnorm/pydamage_results/pydamage_filtered_results.csv", index_col  = 0)

checkm2 = pandas.read_csv("checkm2/quality_report.tsv", sep ="\t", index_col = 0)

positives = set(pydam.reference)

covs = pandas.read_csv("../GenomeBinning/depths/contigs/MEGAHIT-coassembly-bbnorm-depth.txt.gz", compression = "gzip", sep = "\t", index_col=0)
covs = covs[[c for c in covs.columns if not c.endswith("-var")]]
contig_md = covs[["contigLen", "totalAvgDepth"]]
cov_mat = covs[[c for c in covs.columns if c not in contig_md.columns]]
cov_mat.columns = [c.split('-bbnorm-')[1][:-4] for c in cov_mat.columns]
contig_md.loc[:, 'mapped_bases'] = contig_md.contigLen*contig_md.totalAvgDepth
read_mat = cov_mat.apply(lambda x : x*contig_md.contigLen, axis="index")
read_mat.to_csv("contig_read_mat.csv")
cov_mat.to_csv("contig_cov_mat.csv")
real_abs_cont = read_mat.apply(lambda x : x/x.sum(), axis = 0)
relab_ancient = real_abs_cont.loc[positives].sum()
relab_ancient.to_csv("fract_ancient.csv", header = "fraction_ancient")

def bin_stats(seqs) :
    tot_len = sum([len(s) for s in seqs])

    return {
        'size' : tot_len,
        'nb_contigs' : len(seqs),
        'aFreq' : sum([len(s) for s in seqs if s.id in positives])/tot_len,
        'GC' : sum([ s.seq.count("G") + s.seq.count("C") for s in seqs])/tot_len,
        'mapped_bases' :  contig_md.loc[[s.id for s in seqs]].mapped_bases.sum(),
        'mean_cov' : contig_md.loc[[s.id for s in seqs]].mapped_bases.sum() / tot_len ,
        'covs' : read_mat.loc[[s.id for s in seqs]].sum()/tot_len
    }

data = { b[:-3] : bin_stats(list(SeqIO.parse(f"bins/{b}", "fasta"))) for b in os.listdir("bins")}
cov_mat = pandas.DataFrame.from_records({ k : v['covs'] for k,v in data.items()}).transpose()
bin_md = pandas.DataFrame.from_dict({ k : {kk : vv for kk,vv in v.items() if kk != "covs"} for k,v in data.items()}, orient = "index")
bin_md['MAG_quality'] = "LQ"
bin_md.loc[checkm2.index[(checkm2.Completeness > 40) & (checkm2.Contamination < 3)], "MAG_quality"] = 'MQ'
bin_md.loc[checkm2.index[(checkm2.Completeness > 70) & (checkm2.Contamination < 5)], "MAG_quality"] = 'HQ'

cov_mat.to_csv("bin_cov_mat.csv")
bin_md.to_csv("bin_md.csv")
contig_md.to_csv("contig_md.csv")

melted = cov_mat.reset_index().melt(value_vars=cov_mat.columns,  var_name = "sample", value_name = "coverage", id_vars = "index")
melted = pandas.merge(melted, bin_md.reset_index(), on = "index" )

b = 'MEGAHIT-MetaBAT2-coassembly-bbnorm.55'

