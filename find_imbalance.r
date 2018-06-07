suppressMessages(library(plyr))
suppressMessages(library(rmutil))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ape))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))

# check existence of custom function: allele.imbalance(rna, dna, rna.err, dna.err)
source("allele_imbalance.r")

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
dna_file <- args[2]   # whole genome seq counts
rna_file <- args[3]   # exome seq counts
genes_file <- args[4] # gene information, provided by gencode
gq_file <- args[5]    # a tab-delimited file specifying GQ scores for each variant
output_dir <- args[6] # the output directory. must have a trailing slash

import_counts <- function(file_path){
  # import raw counts from WASP
  counts = read.table(gzfile(file_path), sep = " ", header = F , stringsAsFactors= F)
  # name columns such that they can be used in downstream analysis
  colnames(counts) = c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors")
  # add rsID string for later use by allele_imbalance.r
  counts$rsID = paste0(counts$chr, ":", counts$start, "_", counts$ref, "/", counts$alt)
  # remove genotype column. all rows should be heterozygotes already
  counts$genotype = NULL
  counts
}

# import raw dna and rna counts from WASP
dna <- import_counts(dna_file)
rna <- import_counts(rna_file)

# import GQ data from variant_filtration step
gq <- read.table(gq_file, sep="\t", header=T)
# rename GQ column so that it doesn't have the sample name in it
names(gq)[names(gq) == paste0(sample_name, ".GQ")] = "GQ"
# create rsID column for later merging
gq$rsID <- paste0(gq$CHROM, ":", gq$POS, "_", gq$REF, "/", gq$ALT)
gq$genotype.error = 10^{-(gq$GQ/10)}
# retain only the cols that we need
gq = gq[c("genotype.error", "rsID")]

# import gene info from gencode
genes= readGFF(genes_file)


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

g_genes = as(genes, "GRanges")
dna <- proc_counts(dna, g_genes)
rna <- proc_counts(rna, g_genes)
# save space by removing the Granges version of genes. it isn't needed anymore
rm(g_genes)

# add genotype.error from gq to dna by JOINing on common column rsID
dna$genotype.error = gq$genotype.error[match(dna$rsID, gq$rsID)]
# save space by removing gq, it isn't needed anymore
rm(gq)

# intersect dna and rna data frames with each other to get common rows
INT= intersect(dna$rsID, rna$rsID)
dna= dna[dna$rsID %in% INT,]
rna= rna[rna$rsID %in% INT,]

# save the results to csv files
# but make sure the directory exists first!
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = T)
}
write.csv(dna, paste0(output_dir, sample_name, ".dna.csv"), row.names = F)
write.csv(rna, paste0(output_dir, sample_name, ".rna.csv"), row.names = F)

# calculate error rates
err.rate= function(ref, alt, err){
  rate= (sum(err))/sum(ref+alt+err)
  rate*(3/2)
}
dna.err= err.rate(dna$ref.matches, dna$alt.matches, dna$errors)
rna.err= err.rate(rna$ref.matches, rna$alt.matches, rna$errors)

# NOW, we can finally call allele.imbalance()
# note that dna and rna are data frames that need to have these columns:
# 	ref.matches (ie ref_allele_count),
# 	N (total allele count - ref+alt),
# 	genotype.error (ie 10^(-GQ/10) - only appears in dna DF),
# 	rsID (<chr>:<pos>_<allele1>/<allele2>),
# 	gene (name of gene that this SNP appears in, as rnaovided by gencode),
#   start (position of SNP)
dna = dna[c('ref.matches', 'N', 'genotype.error', 'rsID', 'gene', 'start')]
rna = rna[c('ref.matches', 'N', 'rsID', 'gene', 'start')]
# call allele_imbalance and store the result in res
res = as.data.frame(allele.imbalance(rna, dna, rna.err, dna.err))

# rename col 'd' to 'a', since it represents estimates of allelic imbalance
colnames(res)[which(names(res) == "d")] <- "a"
# convert gene col from type factor to type char so dplyr is happy
res$gene = as.character(res$gene)
# add gene names to res. it only has gene_id right now
res = left_join(res, unique(genes[, c("gene_name", "gene_id")]), by = c("gene"="gene_id"))
write.csv(res, paste0(output_dir, sample_name, ".res.csv"), row.names = F)

# # some plots
# hist(res$a, breaks=15, xlab="Estimate of Allelic Imbalance", ylab="Occurence in Genes", main="")
# plot(pval~gene_name, res)

# # plot
# ref= as.numeric(as.character(rna$ref.matches))
# N=rowSums(rna[,c("ref.matches","alt.matches")])
# d= ref/N
# d= d[!is.na(d)]
# plot(density(d), col= "red", main= "density distribution of intersecting sites with ref+alt >= 10", xlab= "ref/N", ylim= c(0,5))
# 
# rd= rna$ref.matches/(rna$N)
# lines(density(rd), col= "blue")
# legend(0.1, 3, legend=c("all", "N>= 10"),
#        col=c("red", "blue"), lty=1, cex=0.8, bty = "n")