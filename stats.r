.libPaths("~/anaconda2/lib/R/library/")
library(plyr)
library(rmutil)
library(GenomicRanges)
library(ape)
library(rtracklayer)
# library(VariantAnnotation)

# custom function: allele.imbalance(rna, dna, rna.err, dna.err)
source("allele_imbalance.r")

setwd("~/allele_specific_analysis/stats")

# sample_name <- args[1]
sample_name <- "Jurkat"

import_counts <- function(sample_name, type){
  # import raw counts from WASP
  counts = read.table(gzfile(paste0(sample_name, ".", type, ".as_counts.txt.gz")), sep = " ", header = F , stringsAsFactors= F)
  # name columns such that they can be used in downstream analysis
  colnames(counts) = c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors")
  # add rsID string for later use by allele_imbalance.r
  counts$rsID = paste0(counts$chr, ":", counts$start, "_", counts$ref, "/", counts$alt)
  # remove genotype column. all rows should be heterozygotes already
  counts$genotype = NULL
  counts
}

# import raw dna and rna counts from WASP
dna <- import_counts(sample_name, "dna")
rna <- import_counts(sample_name, "rna")

# import gene info from gencode
genes= readGFFAsGRanges("gencode.v19.genes.chr.gtf")

# process allele specific counts
proc_counts= function(counts, genes){
  # load counts and plot proportions
  counts$N=rowSums(counts[,c("ref.matches","alt.matches")])
  counts$p= counts$ref.matches/counts$N
  
  # step 1 filter counts
  counts= subset(counts, N >= 10)
  counts$end= counts$start
  counts= GRanges(counts)
  
  # step 2 add genes
  hits= findOverlaps(counts, genes,type = "within")
  names(counts) <- NULL
  counts= as.data.frame(counts[queryHits(hits)])
  genes= as.data.frame(genes[subjectHits(hits)])
  counts$gene= genes$gene_id
  counts= counts[!duplicated(counts),]
}

dna <- proc_counts(dna, genes)
rna <- proc_counts(rna, genes)

# import GQ data from variant_filtration step
gq <- read.table(paste0(sample_name, ".gq_subset.table"), sep="\t", header=T)
# rename GQ column so that it doesn't have the sample name in it
names(gq)[names(gq) == paste0(sample_name, ".GQ")] = "GQ"
# create rsID column for later merging
gq$rsID <- paste0(gq$CHROM, ":", gq$POS, "_", gq$REF, "/", gq$ALT)
gq$genotype.error = 10^{-(gq$GQ/10)}
# retain only the cols that we need
gq = gq[c("genotype.error", "rsID")]

# add genotype.error from gq to dna by JOINing on common column rsID
dna$genotype.error = gq$genotype.error[match(dna$rsID, gq$rsID)]

# intersect dna and rna data frames with each other to get common rows, then save the results to csv files
INT= intersect(dna$rsID, rna$rsID)
dna= dna[dna$rsID %in% INT,]
write.csv(dna, paste0(sample_name, ".dna.csv"), row.names = F)
rna= rna[rna$rsID %in% INT,]
write.csv(dna, paste0(sample_name, ".rna.csv"), row.names = F)

# calculate error rates
err.rate= function(ref, alt, err){
  rate= (sum(err))/sum(ref+alt+err)
  rate*(3/2)
}
dna.err= err.rate(dna$ref.matches, dna$alt.matches, dna$errors)
rna.err= err.rate(rna$ref.matches, rna$alt.matches, rna$errors)

# dna and rna are data frames that need to have these columns:
# 	ref.matches (ie ref_allele_count),
# 	N (total allele count - ref+alt),
# 	genotype.error (ie 10^(-GQ/10) - only appears in dna DF),
# 	rsID (<chr>:<pos>_<allele1>/<allele2>),
# 	gene (name of gene that this SNP appears in, as rnaovided by gencode),
#   start (position of SNP)
dna = dna[c('ref.matches', 'N', 'genotype.error', 'rsID', 'gene', 'start')]
rna = rna[c('ref.matches', 'N', 'rsID', 'gene', 'start')]

res = as.data.frame(allele.imbalance(dna, rna, rna.err, dna.err))
write.csv(res, paste0(sample_name, ".res.csv"), row.names = F)