library(data.table)
library(ggplot2)
library(vegan)
library(pheatmap)
library(ggrepel)

weird_sample = c("sediment_1024")
melted = fread("melted_covs.csv")
contig_cov = fread("contig_cov_mat.csv")
contig_md = fread("contig_md.csv")
bin_md = read.csv("bin_md.csv", row.names=1)
contigs = contig_cov$contigName
contig_cov[,contigName := NULL]
lib_md = fread("~/projects/mosaic/M004_aDNA/data/metadata.csv")

c(sample_1020 = 0.28, sample_1026 = 0.45, sample_1025 = 2.62)

ancient_DNA_freq = read.csv("fract_ancient.csv", row.names=1)
colnames(ancient_DNA_freq) = c('aDNA_fract')
lib_md[, aDNA_freq := ancient_DNA_freq[sample_ngs, "aDNA_fract"] ]
samps = colnames(contig_cov)[!colnames(contig_cov) %in% weird_sample]

sour_dist = dist(t(1-read.csv("~/tmp/aDNA/raws/sourmash_simis.csv")))
sour_mds = metaMDS(sour_dist)
sour_mds = data.table(sour_mds$points)
sour_mds[, sample_ngs := sapply(strsplit(labels(sour_dist),".", fixed = TRUE), head,1)]
sour_mds = sour_mds[lib_md, on = .(sample_ngs)]



bin_covs = read.csv("bin_cov_mat.csv", row.names = 1)

pheatmap(bin_covs, annotation_row=bin_md, annotation_col = ancient_DNA_freq, show_rownames=FALSE, clustering_method = "ward.D2")

mds = metaMDS(t(bin_covs))

mdsP = data.table(mds$points)
mdsP[, sample_ngs := row.names(mds$points)]
mdsP[, aDNA_freq := ancient_DNA_freq[mdsP$sample,]]
mdsP = mdsP[lib_md, on = .(sample_ngs)]


ggplot(sour_mds, aes(x=MDS1, y=MDS2, col=layer, shape=core_id))+geom_point(size = 5)
ggplot(mdsP, aes(x=MDS1, y=MDS2, col=layer, shape=core_id))+geom_point(size = 5)

boreogadus
beluga whale
phytoplantcton
