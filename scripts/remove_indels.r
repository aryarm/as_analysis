# remove indels from a gzipped counts text file and print the result to standard out


# get the input and output files from args
args = commandArgs(trailingOnly = TRUE)
counts = args[1]

# load the hets file from disk (note that it must be gzipped!)
hets = read.table(gzfile(counts), sep=" ", header=F)

# name the cols appropriately
colnames(hets) = c("chr", "start", "ref", "alt", "genotype", "ref.matches", "alt.matches", "errors")

# convert ref and alt cols to character vectors
hets$ref = as.character(hets$ref)
hets$alt = as.character(hets$alt)

# subset
hets = hets[nchar(hets$ref) == 1,]
hets = hets[nchar(hets$alt) == 1,]

# write to stdout
write.table(hets, file=stdout(), quote=FALSE, sep=" ", col.names = F, row.names=F)

# if we got this far, we were probably successful
quit(status=0)