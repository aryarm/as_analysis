# remove indels from a (gzipped!) counts text file and print the result to standard out


# get the input and output files from args
args = commandArgs(trailingOnly = TRUE)
counts = args[1]
separator = args[2]

# load the hets file from disk (note that it must be gzipped!)
hets = read.table(gzfile(counts), sep=separator, header=F)

# name the cols appropriately
colnames(hets)[1:4] = c("chr", "start", "ref", "alt")

# convert ref and alt cols to character vectors
hets$ref = as.character(hets$ref)
hets$alt = as.character(hets$alt)

# subset to keep snps only
hets = hets[nchar(hets$ref) == 1,]
hets = hets[nchar(hets$alt) == 1,]

# write to stdout
write.table(hets, file=stdout(), quote=FALSE, sep=separator, col.names = F, row.names=F)

# if we got this far, we were probably successful
quit(status=0)