library(optparse)
library(futile.logger)
library(stringr)
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-i","--input"),
              help="input file"),
  make_option(c("-o","--output"),
              help = "output"),
  make_option(c("-s","--samplename"),
              help = "sample label")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

##read in csv
flog.info(paste("Reading in input file:",opt$input))
raw.csv <- read.csv(opt$input,sep="\t",header=F,quote="")
unmpliup <- raw.csv[raw.csv[,4] == 0,]
rownums <- nrow( unmpliup )
unmpliup[,5] <- rep("unknow",rownums)
unmpliup[,6] <- rep("unknow",rownums)
analysis <- raw.csv[ raw.csv[,4] != 0, ]

analysis <- rbind( analysis, unmpliup )
colnames(analysis) <- analysis[1,]
analysis <- analysis[-1,]

analysis <- analysis[analysis[,1] !="",]



write.table(analysis, opt$output, row.names=F, sep="\t", quote=F)
