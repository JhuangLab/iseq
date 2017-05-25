library(optparse)
library(futile.logger)
library(stringr)
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-i","--input"),
              help="input fmt result file"),
  make_option(c("-a","--casesnv"),
              help="input case snv frq file"),
  make_option(c("-b","--caseindel"),
              help="input case indel frq file"),
  make_option(c("-c","--controlsnv"), default="",
              help="input control snv frq file"),
  make_option(c("-d","--controlindel"),default="",
              help="input control indel frq file"),
  make_option(c("-o","--output"),
              help = "output"),
  make_option(c("-s","--samplename"),
              help = "sample label")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

##read in csv
result <- read.csv(opt$input,sep="\t",header=T,quote = "")
isSNV <- !(str_detect(result$Base.change,"-")) 
#case data 
casesnvFrq <- read.csv(opt$casesnv,sep="\t",header=T,quote = "")
caseindelFrq <- read.csv(opt$caseindel,sep="\t",header=T,quote = "")
chrpos_result <- paste0( result$Chr,":", result$Start )
chrpos_snv <- paste0( casesnvFrq$Chr,":", casesnvFrq$Postion )
chrpos_indel <- paste0( caseindelFrq$Chr,":", caseindelFrq$Postion )
#case snv
snvindex <- match( chrpos_result[isSNV], chrpos_snv )
result[isSNV,"caseAlleDepth"] <- casesnvFrq[snvindex,"AlleDepth"]
result[isSNV,"caseCommaNum"] <- casesnvFrq[snvindex,"CommaNum"]
result[isSNV,"caseDotNum"] <- casesnvFrq[snvindex,"DotNum"]
result[isSNV,"caseMutationType"] <- casesnvFrq[snvindex,"Alt"]
result[isSNV,"caseMutationDepth"] <- casesnvFrq[snvindex,"AltNum"]
result[isSNV,"caseMutationFrq"] <- casesnvFrq[snvindex,"Frquency"]
#case indel
indelindex <- match( chrpos_result[!isSNV], chrpos_indel )
result[!isSNV,"caseAlleDepth"] <- caseindelFrq[indelindex,"AlleDepth"]
result[!isSNV,"caseCommaNum"] <- caseindelFrq[indelindex,"CommaNum"]
result[!isSNV,"caseDotNum"] <- caseindelFrq[indelindex,"DotNum"]
result[!isSNV,"caseMutationType"] <- caseindelFrq[indelindex,"Alt"]
result[!isSNV,"caseMutationDepth"] <- caseindelFrq[indelindex,"AltNum"]
result[!isSNV,"caseMutationFrq"] <- caseindelFrq[indelindex,"Frquency"]

if (opt$controlsnv != ""){
    #control data
    controlsnvFrq <- read.csv(opt$controlsnv,sep="\t",header=T,quote = "")
    controlindelFrq <- read.csv(opt$controlindel,sep="\t",header=T,quote = "")
    chrpos_result <- paste0( result$Chr,":", result$Start )
    chrpos_snv <- paste0( controlsnvFrq$Chr,":", controlsnvFrq$Postion )
    chrpos_indel <- paste0( controlindelFrq$Chr,":", controlindelFrq$Postion )

    #control snv
    snvindex <- match( chrpos_result[isSNV], chrpos_snv )
    result[isSNV,"controlAlleDepth"] <- controlsnvFrq[snvindex,"AlleDepth"]
    result[isSNV,"controlCommaNum"] <- controlsnvFrq[snvindex,"CommaNum"]
    result[isSNV,"controlDotNum"] <- controlsnvFrq[snvindex,"DotNum"]
    result[isSNV,"controlMutationType"] <- controlsnvFrq[snvindex,"Alt"]
    result[isSNV,"controlMutationDepth"] <- controlsnvFrq[snvindex,"AltNum"]
    result[isSNV,"controlMutationFrq"] <- controlsnvFrq[snvindex,"Frquency"]
    #control indel
    indelindex <- match( chrpos_result[!isSNV], chrpos_indel )
    result[!isSNV,"controlAlleDepth"] <- controlindelFrq[indelindex,"AlleDepth"]
    result[!isSNV,"controlCommaNum"] <- controlindelFrq[indelindex,"CommaNum"]
    result[!isSNV,"controlDotNum"] <- controlindelFrq[indelindex,"DotNum"]
    result[!isSNV,"controlMutationType"] <- controlindelFrq[indelindex,"Alt"]
    result[!isSNV,"controlMutationDepth"] <- controlindelFrq[indelindex,"AltNum"]
    result[!isSNV,"controlMutationFrq"] <- controlindelFrq[indelindex,"Frquency"]
}
isdel <- (!is.na(result$Mutation.type)) & str_detect(result$Base.change,">-")
result[isdel,"Start"] <- as.numeric(result[isdel,"Start"]) + 1
write.table(result,opt$output,row.names=F,sep="\t",quote=F)

fil <- str_detect(result$Func.refGene,"splicing|exonic") & (is.na(result$Mutation.type)|result$Mutation.type == "NA" | str_detect(result$Mutation.type,"nonsynonymous|frame|stop|missense|nonsense|unknow"))
exon.out <- result[fil,]

out.txt = as.character(opt$output)
exon.out.txt = str_replace(out.txt,'.txt','.exon.txt')
write.table(exon.out,exon.out.txt,sep="\t",col.names=T,row.names=F,quote=F)
