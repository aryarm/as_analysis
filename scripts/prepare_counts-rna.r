suppressMessages(library(plyr))
suppressMessages(library(rmutil))
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressMessages(library(dplyr))
suppressMessages(library(tools))

args = commandArgs(trailingOnly = TRUE)
# import counts for this sample
counts = read.table(gzfile(args[1]), sep=" ", header=F, stringsAsFactors=F, col.names=c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors"))
# import gq data for this sample
use_default_gq = !is.na(as.integer(args[2]))
if (use_default_gq) {
  gq = as.integer(args[2])
} else {
  gq = read.table(gzfile(args[2]), sep="\t", header=F, col.names=c("CHROM", "POS", "REF", "ALT", "GQ"))
}
# import gene info from gencode
if (file_ext(args[3]) == "gtf" | file_ext(args[3]) == "gff") {
  targets = as(readGFF(args[3]), "GRanges")
} else if (file_ext(args[3]) == "bed") {
  targets = import(args[3], format = "bed")
} else if (file_ext(args[3]) == "narrowPeak") {
  targets = import(args[3], format = "narrowPeak")
} else {
  stop("Aborting! Targets file extension is not supported. Must be one of gtf, gff, bed, or narrowPeak.")
}

# what is the output directory prefix? note that it should have a trailing slash
output_dir = args[4]

# switch ref and alt counts if genotype is 1|0
counts = transform(counts, ref.matches = ifelse(genotype == "1|0", alt.matches, ref.matches), alt.matches = ifelse(genotype == "1|0", ref.matches, alt.matches))
# create unique rsID column for later merging with gq
counts$rsID = paste0(counts$chr, ":", counts$start, "_", counts$ref, "/", counts$alt)
# remove genotype col, since all rows should be heterozygotes already
counts$genotype = NULL

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
proc_counts= function(counts, targets){
  # load counts and plot proportions
  counts$N=rowSums(counts[,c("ref.matches","alt.matches")])
  counts$ref_prop= counts$ref.matches/counts$N
  
  # step 1 filter counts
  message("- Removing ", sum(counts$N < 10), " SNPs that have low read counts")
  counts= subset(counts, N >= 10)
  num_old_counts = nrow(counts)
  if (num_old_counts == 0) {
    stop("Aborting! There aren't any SNPs with a sufficient number of supporting reads.")
  }
  counts$end= counts$start
  counts= GRanges(counts)

  message("Make sure these match up!\nTargets seqStyle: ", seqlevelsStyle(targets), "\nCounts seqStyle: ", seqlevelsStyle(counts))
  
  # step 2 add targets
  hits= findOverlaps(counts, targets, type = "within")
  names(counts) = NULL
  counts= as.data.frame(counts[queryHits(hits)])
  targets= as.data.frame(targets[subjectHits(hits)])
  if ("gene_id" %in% colnames(targets)){
    counts$target = targets$gene_id
  } else {
    counts$target = paste(paste(targets$seqnames, targets$start, sep=":"), targets$end, sep="-")
  }
  counts= counts[!duplicated(counts),]
  message("- Removed ", num_old_counts-nrow(counts), " SNPs that don't overlap a target region")
  if (nrow(counts) == 0) {
    stop("Aborting! There aren't any SNPs that lie within the target regions.")
  }
  counts
}

# process count data, adding gene information
message("Processing read counts of ", nrow(counts), " SNPs...")
counts = proc_counts(counts, targets)
# save memory by removing genes. it isn't needed anymore
rm(targets)

# add genotype.error from gq to rna by JOINing on common column rsID
message("Adding genotype error to ", nrow(counts), " RNA SNPs...")
if (use_default_gq) {
  counts$genotype.error = gq
} else {
  counts = inner_join(counts, gq, by="rsID")
  # save memory by removing gq, it isn't needed anymore
  rm(gq)
}

# save the results to csv files
message("Writing ", nrow(counts), " final SNPs to file...\n")
counts_output = paste0(output_dir, "rna.csv")
write.csv(counts, counts_output, row.names = F, quote = F)
system(paste("gzip", counts_output))

# if we got this far, we were probably successful
quit(status=0)
