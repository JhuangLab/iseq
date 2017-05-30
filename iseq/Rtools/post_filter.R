library(optparse)
library(futile.logger)
library(stringr)
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-i","--inputcsvs"),
              help="input csv"),
  make_option(c("--min-casedepth"),dest="mincasedepth",
              help="Mininm case total depth",default=0),
  make_option(c("--min-casefrq"),dest="mincasefrq",
              help="Mininm case mutation frquency",default=0.00),
  make_option(c("--max-controlfrq"),dest="maxcontrolfrq",
              help="Maxinm control mutation frquency",default=1),
  make_option(c("--mapper_field"),dest="mapper_field",default=3,
              help="File name of which field is mapper name , eg. L01.exon.mapper.caller.txt, mapper field is 3"),
  make_option(c("--caller_field"),dest="caller_field",default=4,
              help="File name of which field is caller name , eg. L01.exon.mapper.caller.txt, caller field is 4"),
  make_option(c("-o","--outputcsv"),
              help = "output csv"),
  make_option(c("-s","--samplename"),
              help = "sample label"),
  make_option(c("-m","--mode"),default="somatic",
              help = "somatic or germline")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
files <- str_split(opt$inputcsvs,",")[[1]]
flog.info(paste("Reading in csv files:",opt$inputcsv))
csv.out <- NULL
for(file in files){
    if(file.exists(file)){
        raw.csv <- read.table(file,sep="\t",header=T)
    }
    else{
        next
    }
    mapper <- str_split(file,fixed("."))[[1]][opt$mapper_field]
    caller <- str_split(file,fixed("."))[[1]][opt$caller_field]
    mapper <- rep(mapper,nrow(raw.csv)) 
    caller <- rep(caller,nrow(raw.csv)) 
    raw.csv <- cbind(raw.csv,mapper,caller)
    csv.out <- rbind(csv.out,raw.csv)
}
filcasedepth <- as.numeric(csv.out$caseAlleDepth) >= opt$mincasedepth
filcasefrq <- as.numeric(csv.out$caseMutationFrq) >= opt$mincasefrq
if( opt$mode == "somatic" ){
    filcontrolfrq <- as.numeric(csv.out$controlMutationFrq) <= opt$maxcontrolfrq
    issame <- str_detect(fixed(csv.out$controlMutationType), fixed(csv.out$caseMutationType))
    postfil <- (filcasedepth & filcasefrq & (filcontrolfrq | !issame )) | (csv.out$Gene == "FLT3" && str_detect(csv.out$Mutation.type,"frame"))
}else{
    postfil <- (filcasedepth & filcasefrq)
} 
getmapcaller <- function(csv.out){
    chrpos <- paste0(csv.out$Patient.ID,":",csv.out$Chr,":",csv.out$Start)
    chrpos <- unique(chrpos)
    result <- NULL
    for(i in chrpos){
        fil <- paste0(csv.out$Patient.ID,":",csv.out$Chr,":",csv.out$Start) == i
        tmp <- csv.out[fil,] 
        mapper = NULL
        for(j in unique(tmp[,"mapper"])){
            if( is.null(mapper) ){
                mapper <- j
            }else{
                mapper <- paste0(mapper,",",j)
            }
        }
        caller = NULL
        for(j in unique(tmp[,"caller"])){
            if(is.null(caller)){
                caller <- j
            }else{
                caller <- paste0(caller,",",j)
            }
        }
        tmp[,"mapper"] <- mapper
        tmp[,"caller"] <- caller
        result <- rbind(result,tmp[1,])
    }
    result <- result[!is.na(result[,1]),]
    return(result)
}
sortcaller <- function(caller){
    result <- ""
    order <- "Pindel,Tvc,Varscan,Mutect,Lofreq,Haplotypecaller,Unifiedgenotyper"
    order <- str_split(order,",")[[1]]
    for(i in order){
       if(str_detect(caller,i)){
           if(result == ""){
               result <- i
           }else{
               result <- paste0(result,",",i)
           }
       } 
    }
    return(result)
}

fil2anaysisIndel <- (is.na(csv.out$caseMutationFrq) | (csv.out$caseMutationFrq == 0) ) & (str_detect(csv.out$Base.change,"-")|str_length(csv.out$Base.change)>3)
csv.out.indel <- csv.out[fil2anaysisIndel,]
if(nrow(csv.out.indel)>0){
    csv.out.indel <- getmapcaller( csv.out.indel )
    for(i in 1:nrow(csv.out.indel)){
        csv.out.indel[i,"caller"] <- sortcaller(csv.out.indel[i,"caller"])
    }
}
flog.info(paste("Writing in csv files:",paste0(as.character(opt$outputcsv),".del")))
write.table(csv.out.indel, paste0(as.character(opt$outputcsv),".del"),sep="\t",col.names=T,row.names=F,quote=F)
csv.out <- csv.out[postfil,]
csv.out <- getmapcaller( csv.out ) 
if(nrow(csv.out) == 0){
  q()
}
for(i in 1:nrow(csv.out)){
    csv.out[i,"caller"] <- sortcaller(csv.out[i,"caller"])
}
flog.info(paste("Writing in csv files:",opt$outputcsv))
write.table(csv.out, as.character(opt$outputcsv),sep="\t",col.names=T,row.names=F,quote=F)
