suppressMessages(library(plyr))
suppressMessages(library(rmutil))
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly = TRUE)
# import dna counts for this sample
dna = read.table(gzfile(args[1]), sep=" ", header=F, stringsAsFactors=F, col.names=c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors"))
# import rna counts for this sample
rna = read.table(gzfile(args[2]), sep=" ", header=F, stringsAsFactors=F, col.names=c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors"))
# import gq data for this sample
use_default_gq = !is.na(as.integer(args[3]))
if (use_default_gq) {
  gq = as.integer(args[3])
} else {
  gq = read.table(gzfile(args[3]), sep="\t", header=F, col.names=c("CHROM", "POS", "REF", "ALT", "GQ"), stringsAsFactors=F)
}
# import gene info from gencode
genes= as(readGFF(args[4]), "GRanges")
# what is the output directory prefix? note that it should have a trailing slash
output_dir = args[5]

# switch ref and alt counts if genotype is 1|0
dna = transform(dna, ref.matches = ifelse(genotype == "1|0", alt.matches, ref.matches), alt.matches = ifelse(genotype == "1|0", ref.matches, alt.matches))
rna = transform(rna, ref.matches = ifelse(genotype == "1|0", alt.matches, ref.matches), alt.matches = ifelse(genotype == "1|0", ref.matches, alt.matches))
# create unique rsID column for later merging with gq
dna$rsID = paste0(dna$chr, ":", dna$start, "_", dna$ref, "/", dna$alt)
rna$rsID = paste0(rna$chr, ":", rna$start, "_", rna$ref, "/", rna$alt)
# remove genotype col, since all rows should be heterozygotes already
dna$genotype = NULL
rna$genotype = NULL

# process gq data
if (use_default_gq) {
  # calculate the genotype error
  gq = 10^{-(gq/10)}
} else {
  # create rsID column for later merging and get rid of other columns
  gq$rsID = paste0(gq$CHROM, ":", gq$POS, "_", gq$REF, "/", gq$ALT)
  gq$CHROM = NULL
  gq$POS = NULL
  gq$REF = NULL
  gq$ALT = NULL
  # remove any SNPs for which the genotype quality is unknown
  # and then convert it to a numeric type
  gq = gq[gq$GQ != ".",]
  gq$GQ = as.numeric(as.character(gq$GQ))
  # generate the genotype.error column
  gq$genotype.error = 10^{-(gq$GQ/10)}
  # retain only the cols that we need
  gq = gq[c("genotype.error", "rsID")]
}

# process allele specific counts
proc_counts= function(counts, genes){
  # load counts and proportions
  counts$N=rowSums(counts[,c("ref.matches","alt.matches")])
  counts$p= counts$ref.matches/counts$N
  
  # step 1 filter counts
  message("- Removing ", sum(counts$N < 10), " SNPs that have low read counts")
  counts= subset(counts, N >= 10)
  num_old_counts = nrow(counts)
  counts$end= counts$start
  counts= GRanges(counts)
  
  # step 2 add genes
  hits= findOverlaps(counts, genes, type = "within")
  names(counts) = NULL
  counts= as.data.frame(counts[queryHits(hits)])
  genes= as.data.frame(genes[subjectHits(hits)])
  counts$gene= genes$gene_id
  counts= counts[!duplicated(counts),]
  message("- Removed ", num_old_counts-nrow(counts), " SNPs that have don't overlap a gene")
  counts
}

# process dna and rna count data, adding gene information
message("Processing dna counts of ", nrow(dna), " SNPs...")
dna = proc_counts(dna, genes)
message("Processing rna counts of ", nrow(rna), " SNPs...")
rna = proc_counts(rna, genes)
# save memory by removing genes. it isn't needed anymore
rm(genes)

# add genotype.error from gq to dna by JOINing on common column rsID
message("Adding genotype error to ", nrow(dna), " DNA SNPs...")
if (use_default_gq) {
  dna$genotype.error = gq
} else {
  dna = inner_join(dna, gq, by="rsID")
  # save memory by removing gq, it isn't needed anymore
  rm(gq)
}

# intersect dna and rna data frames with each other to get common rows
message("Removing any SNPS that aren't present in both DNA and RNA because of filtering...")
INT= intersect(dna$rsID, rna$rsID)
dna= dna[dna$rsID %in% INT,]
rna= rna[rna$rsID %in% INT,]

# save the results to csv files
message("Writing ", nrow(dna), " final SNPs to file...\n")
dna_output = paste0(output_dir, "dna.csv")
rna_output = paste0(output_dir, "rna.csv")
write.csv(dna, dna_output, row.names = F)
system(paste("gzip", dna_output))
write.csv(rna, rna_output, row.names = F)
system(paste("gzip", rna_output))

# if we got this far, we were probably successful
quit(status=0)
