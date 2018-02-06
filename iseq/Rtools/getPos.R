library(optparse)
library(futile.logger)
library(stringr)
library(data.table)
options(stringsAsFactors = FALSE)
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE, 
  help = "Print extra output [default]"), make_option(c("-i", "--input"), help = "input file"), 
  make_option(c("-o", "--output"), help = "output"), make_option(c("-s", "--samplename"), 
    help = "sample label"), make_option(c("-e", "--exononly"), action = "store_true", 
    default = FALSE, help = "only get exon pos"), make_option(c("-p", "--split"), 
    help = "split character"))
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list))

## read in csv
flog.info(paste("Reading in input file:", opt$input))
raw.csv <- fread(opt$input, sep = opt$split, header = T)
raw.csv <- as.data.frame(raw.csv)
if (opt$exononly) {
  # fil <- str_detect(raw.csv$Func.refGene,'exon|splic|stream|stop')
  fil <- !str_detect(raw.csv$Func.refGene, "intron|inte")
  raw.csv <- raw.csv[fil, ]
}
isdel <- (!is.na(raw.csv$ExonicFunc.refGene)) & str_detect(raw.csv$Alt, "-")
raw.csv[isdel, "Start"] <- as.numeric(raw.csv[isdel, "Start"]) - 1
pos <- cbind(raw.csv$Chr, raw.csv$Start)
pos <- pos[pos[, 1] != "contig", ]
colnames(pos) = pos[1, ]
# pos = pos[-1,]
fwrite(pos, opt$output, row.names = F, sep = "\t", quote = F)
