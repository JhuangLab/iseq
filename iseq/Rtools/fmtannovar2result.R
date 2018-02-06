library(optparse)
library(futile.logger)
library(stringr)
library(data.table)
library(ngstk)
options(stringsAsFactors = FALSE)
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE, 
  help = "Print extra output [default]"), make_option(c("-i", "--inputcsv"), help = "input csv"), 
  make_option(c("-o", "--outputcsv"), help = "output csv"), make_option(c("-s", 
    "--samplename"), help = "sample label"), make_option(c("-c", "--colnames"), 
    help = "colnames"))
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list))
tb.colnames <- colnames(fread(sprintf("head -n 1 %s", opt$inputcsv)))
handler_fun <- function(x = "", i = 1, ...) {
  ## read in csv
  #flog.info(paste("Reading in csv file:", opt$inputcsv))
  #raw.csv <- fread(opt$inputcsv, sep = ",", header = T)
  #raw.csv <- as.data.frame(raw.csv)
  raw.csv <- as.data.frame(x)
  raw.csv <- raw.csv[raw.csv[, 1] != "contig", ]
  if (i != 1) {
    colnames(raw.csv) <- tb.colnames
  } else {
    raw.csv <- raw.csv[-1,]
    colnames(raw.csv) <- tb.colnames
  }
  colnames.str = opt$colnames
  colnames.need <- str_split(colnames.str, ",")[[1]]
  
  print(head(raw.csv$Chr))
  fil <- !str_detect(raw.csv$Chr, "chrM|chrUn|random")
  csv.out <- raw.csv[fil, ]
  csv.out <- raw.csv[, colnames(csv.out) %in% colnames.need]
  rm(raw.csv)
  gc()
  csv.out[csv.out == "."] = NA
  ## output
  out.total = nrow(csv.out)
  flog.info(paste(out.total, "mutations are output to file:", opt$outputcsv))
  Patient.ID = rep(as.character(opt$samplename), out.total)
  Base.change = paste(as.character(csv.out$Ref), as.character(csv.out$Alt), sep = ">")
  print(opt$samplename)
  print(out.total)
  csv.out <- data.frame(Patient.ID, csv.out)
  Ref.col <- which(colnames(csv.out) == "Ref")
  Alt.col <- which(colnames(csv.out) == "Alt")
  csv.out <- data.frame(csv.out[1:(Ref.col - 1)], Base.change, csv.out[, (Alt.col + 
    1):ncol(csv.out)])
  index <- is.null(csv.out$AAChange.refGene) | is.na(csv.out$AAChange.refGene) & !is.na(csv.out$GeneDetail.refGene)
  csv.out$AAChange.refGene[index] <- csv.out$GeneDetail.refGene[index]
  csv.out <- csv.out[, -(which(colnames(csv.out) == "GeneDetail.refGene"))]
  colnames(csv.out)[which(colnames(csv.out) == "Gene.refGene")] <- "Gene"
  colnames(csv.out)[which(colnames(csv.out) == "ExonicFunc.refGene")] <- "Mutation.type"
  colnames(csv.out)[which(colnames(csv.out) == "SIFT_pred")] <- "SIFT.prediction"
  colnames(csv.out)[which(colnames(csv.out) == "SIFT_score")] <- "SIFT.score"
  colnames(csv.out)[which(colnames(csv.out) == "Polyphen2_HDIV_score")] <- "PolyPhen.score"
  colnames(csv.out)[which(colnames(csv.out) == "Polyphen2_HDIV_pred")] <- "PolyPhen.prediction"
  colnames(csv.out)[which(colnames(csv.out) == "AAChange.refGene")] <- "Amino.acid.change"
  replace_gene_fun <- function(x, y) {
    return(str_replace_all(x, paste0(y, ":"), ""))
  }
  flog.info("Processing Amino.acid.change...")
  csv.out$Amino.acid.change <- apply(csv.out, 1, function(x) {
    replace_gene_fun(x["Amino.acid.change"], 
      x["Gene"])
  })
  flog.info("Amino.acid.change processing finised.")
  ## change Mutation type
  csv.out$Mutation.type = str_replace_all(csv.out$Mutation.type, "nonsynonymous SNV", 
    "missense")
  csv.out$Mutation.type = str_replace_all(csv.out$Mutation.type, "stopgain", "nonsense")
  csv.out$Mutation.type = str_replace_all(csv.out$Mutation.type, "insertion", "ins")
  csv.out$Mutation.type = str_replace_all(csv.out$Mutation.type, "deletion", "del")
  csv.out$SIFT.prediction = str_replace(csv.out$SIFT.prediction, "T", "Tolerated")
  csv.out$SIFT.prediction = str_replace(csv.out$SIFT.prediction, "D", "Damaging")
  csv.out$PolyPhen.prediction = str_replace(csv.out$PolyPhen.prediction, "B", "benign")
  csv.out$PolyPhen.prediction = str_replace(csv.out$PolyPhen.prediction, "D", "probably damaging")
  csv.out$PolyPhen.prediction = str_replace(csv.out$PolyPhen.prediction, "P", "possibly damaging")
  csv.out.bak <- csv.out
  isdel <- (!is.na(csv.out$Mutation.type)) & str_detect(csv.out$Base.change, ">-")
  csv.out[isdel, "Start"] <- as.numeric(csv.out[isdel, "Start"]) - 1
  if (i == 1) {
    fwrite(csv.out, as.character(opt$outputcsv), sep = "\t", col.names = T, row.names = F, 
      quote = F)
  } else {
    fwrite(csv.out, as.character(opt$outputcsv), sep = "\t", col.names = F, row.names = F, 
      quote = F, append = TRUE)
  }
}

batch_file(opt$inputcsv, 500003, handler_fun, extra_params = list(tb.colnames = tb.colnames, opt = opt), 
           extra_fread_params = list(sep = ",", header = FALSE, return_1L = FALSE))
