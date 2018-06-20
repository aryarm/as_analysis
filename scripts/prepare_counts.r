suppressMessages(library(plyr))
suppressMessages(library(rmutil))
suppressMessages(library(ape))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly = TRUE)
# The path to a text file specifying where to find dna/rna count and gq files for each sample.
# Each row in the sample file should represent a different sample.
# The sample file should have 4 columns (each separated by a single tab):
#       <unique_sample_name> <dna_counts_path> <rna_counts_path> <genotype_qual_path>
samples = read.table(args[1], sep="\t", header=F)
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