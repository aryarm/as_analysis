suppressMessages(library(plyr))
suppressMessages(library(rmutil))
suppressMessages(library(ape))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))

# check existence of custom function: allele.imbalance(rna, dna, rna.err, dna.err)
source("allele_imbalance.r")

args = commandArgs(trailingOnly = TRUE)
# The path to a text file specifying where to find dna/rna count and gq files for each sample.
# Each row in the sample file should represent a different sample.
# The sample file should have 4 columns (each separated by a single space):
#       <unique_sample_name> <dna_counts_path> <rna_counts_path> <genotype_qual_path>
samples = read.table(args[1], sep=" ", header=F)
colnames(samples) = c("name", "dna", "rna", "gq")
# import gene info from gencode
genes= readGFF(args[2])
g_genes = as(genes, "GRanges")
# what is the output directory prefix? note that it must have a trailing slash
output_dir = args[3]

import_counts = function(sample_paths){
  # import raw counts from WASP
  counts = ldply(apply(sample_paths, MARGIN=1, function(sample){
    counts = read.table(gzfile(sample[2]), sep=" ", header=F, stringsAsFactors=F)
    counts$sample = rep(sample[1], nrow(counts))
    counts
  }), data.frame)
  # name columns such that they can be used in downstream analysis
  colnames(counts) = c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors", "sample")
  # add rsID string for later use by allele_imbalance.r
  counts$rsID = paste0(counts$chr, ":", counts$start, "_", counts$ref, "/", counts$alt)
  # remove genotype col, since all rows should be heterozygotes already
  counts$genotype = NULL
  counts
}

dna = import_counts(samples[c("name", "dna")])
rna = import_counts(samples[c("name", "rna")])

# import raw qg data from GATK
gq = ldply(apply(samples[c("name", "gq")], MARGIN=1, function(sample){
  gq = read.table(sample[2], sep="\t", header=T)
  # rename GQ column so that it doesn't have the sample name in it
  names(gq)[names(gq) == paste0(sample[1], ".GQ")] = "GQ"
  gq$sample = rep(sample[1], nrow(gq))
  gq
}), data.frame)
# create rsID column for later merging
gq$rsID = paste0(gq$CHROM, ":", gq$POS, "_", gq$REF, "/", gq$ALT)
gq$genotype.error = 10^{-(gq$GQ/10)}
# retain only the cols that we need
gq = gq[c("genotype.error", "rsID", "sample")]


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
  hits= findOverlaps(counts, genes, type = "within")
  names(counts) = NULL
  counts= as.data.frame(counts[queryHits(hits)])
  genes= as.data.frame(genes[subjectHits(hits)])
  counts$gene= genes$gene_id
  counts= counts[!duplicated(counts),]
}

# process dna and rna count data, adding gene information
dna = proc_counts(dna, g_genes)
rna = proc_counts(rna, g_genes)
# save space by removing the Granges version of genes. it isn't needed anymore
rm(g_genes)

# add genotype.error from gq to dna by JOINing on common column rsID
dna = merge(dna, gq, by=c("sample", "rsID"), sort=F)
# save space by removing gq, it isn't needed anymore
rm(gq)

# the rest of the pipeline is applied per sample
per_sample = function(sample_name, dna, rna) {
  # change output_dir for this sample
  output_dir = paste0(output_dir, sample_name)
  
  # intersect dna and rna data frames with each other to get common rows
  INT= intersect(dna$rsID, rna$rsID)
  dna= dna[dna$rsID %in% INT,]
  rna= rna[rna$rsID %in% INT,]
  
  # save the results to csv files
  # but make sure the directory exists first!
  if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }
  write.csv(dna, paste0(output_dir, "/dna.csv"), row.names = F)
  write.csv(rna, paste0(output_dir, "/rna.csv"), row.names = F)
  
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
  colnames(res)[which(names(res) == "d")] = "a"
  # convert gene col from type factor to type char so dplyr is happy
  res$gene = as.character(res$gene)
  # add gene names to res. it only has gene_id right now
  res = left_join(res, unique(genes[, c("gene_name", "gene_id")]), by = c("gene"="gene_id"))
  write.csv(res, paste0(output_dir, "/result.csv"), row.names = F)
  
  # # some plots
  # hist(res$a, breaks=15, xlab="Estimate of Allelic Imbalance", ylab="Occurence in Genes", main="")
  # plot(pval~gene_name, res)
}

# get list of sample names
sample_names= unique(dna$sample)
# split dna and rna by sample, creating a list of dna/rna data frames for each sample
dna_samples = split(dna, f = factor(dna$sample, levels=sample_names))
rna_samples = split(rna, f = factor(rna$sample, levels=sample_names))
# complete the per-sample pipeline with each dna/rna pair of data frames
invisible(mapply(per_sample, sample_name=sample_names, dna=dna_samples, rna=rna_samples))

# if we got this far, we were probably successful
quit(status=0)