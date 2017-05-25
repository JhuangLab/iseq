#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2015 Jinyan HUANG <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License.
@RNA-seq pipeline
@status:  experimental
@version: 1.0 
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from preprocess import *
from variantcaller import *
from refinement import *
# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level = 20,
                    format = '%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt = '%a, %d %b %Y %H:%M:%S',
                    stream = sys.stderr,
                    filemode = "w"
                   )

# ------------------------------------
# Misc functions
# ------------------------------------
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info
# ------------------------------------
# Debug file
# ------------------------------------
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = "usage: %prog -c config.cfg -s A01A -m Fastq2vcf -1 A01_1.fq -2 A02_2.fq --bamprocess 00101111 -o vcf"
    description = "Please set the sample name. e.g. L04A, L04C, L04T."
    optparser = OptionParser(version = "rnaseq 1.0", description = description, usage = usage, add_help_option = False)
    optparser.add_option("-h", "--help", action = "help", help = "Show this help message and exit.")
    optparser.add_option("-c", "--config", dest = "config", default = "config.cfg" ,type = "string",
                         help = "Set the config File.[config.cfg]")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.(Required)")
    optparser.add_option("-m", "--mode", dest = "mode" ,type = "string",
                         help = "Run mode, [genome_index, Fastq2vcf, Fastq2bam, Bam2vcf, Bamprocess].(Required)")
    optparser.add_option("-1", "--fastq1", dest = "fastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-2", "--fastq2", dest = "fastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-d", "--bamprocess", dest = "bamprocess", default = "00101111",type = "string",
                         help = "Need point 8 digit number, eg. 00101111:  Add Read Group, Reorder Contig, Mark Duplicates, SplitNtrim ,RealignerTargetCreator,IndelRealigner,Recalibration,PrintReads, step one by one;01000000 only conduct Reorder contig; 00010000:Only conduct SplitNtrim step.[00101111].If --mode=Bamprocess,this parameter is required")
    optparser.add_option("-i", "--in_bam", dest = "in_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 and is empty, the in_bam is required!)")
    optparser.add_option("-o", "--out_dir", dest = "out_dir" ,type = "string", default = "vcf",
                         help = "Set the vcf file out_dir.[vcf]")
    return(optparser)
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options, args) = optparser.parse_args()
    if not options.config:
        optparser.print_help()
        sys.exit(1)
    elif not options.samplename:
        optparser.print_help()
        sys.exit(1)
    elif (not options.fastq1 and not options.fastq2 and not options.in_bam) and options.mode != "Genomeindex":
        optparser.print_help()
        sys.exit(1)
    elif not options.mode:  
        optparser.print_help()
        sys.exit(1)
    elif options.mode.lower() not in ["genomeindex", "fastq2vcf", "fastq2bam", "bam2vcf", "bamprocess"]:  
        optparser.print_help()
        sys.exit(1)
    elif options.mode in ["Fastq2vcf","Fastq2bam"] and not options.fastq1 and not options.fastq2: 
        optparser.print_help()
        sys.exit(1)
    elif options.mode in ["Bam2vcf","Bamprocess"] and not options.in_bam: 
        optparser.print_help()
        sys.exit(1)
    return(options)

def panel(options=""):
    if options == "":
        options = opt_validate(prepare_optparser())
    cfg = get_config(options.config)
    vcfout_dir = options.out_dir
    ################ Genome Index Mode #################
    options.mode = options.mode.lower()
    if options.mode == "genomeindex":
        options.genome_index = "1"
        options.fastq_mapping = "0"
        options.bamprocess = "00000000"
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"

    ################ Fastq2vcf Mode ################### 
    if options.mode == "fastq2vcf":
        options.genome_index = "0"
        options.fastq_mapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = ""
        options.vcfannovar = "1"
        options.mpileup = "1"
    ################ Fastq2bam Mode ################### 
    if options.mode == "fastq2bam":
        options.genome_index = "0"
        options.fastq_mapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"
        options.mpileup = "0"
    ################ Bam2vcf Mode ################### 
    if options.mode == "bam2vcf":
        options.genome_index = "0"
        options.fastq_mapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "" 
        options.vcfannovar = "1"
        options.mpileup = "1"
    if options.mode == "bamprocess":
        options.genome_index = "0"
        options.fastq_mapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "" 
        options.vcfannovar = "0"
        options.mpileup = "0"

    ################ Main Region ###################### 
    bamfile_list = pre_process(options)
    options.seq_type = "dna"
    if options.variantcaller == "1":
        for key in bamfile_list.keys():
            bamfile = bamfile_list[key]
            options.in_bam = bamfile
            options.out_dir = vcfout_dir + "/" + key
            vcffile_list = variant_caller(options)
            for vcf in vcffile_list.keys():
                key = key.capitalize()
                vcffile = vcffile_list[vcf]
                if vcffile:
                    options.out_dir = vcffile.dirname + "/" 
                    options.case_vcf = str(vcffile)
                    vcf = vcf.capitalize()
                    if vcf in "UnifiedGenotyper":
                        options.vcffilter = "wespipeline"
                    options.exononly = True
                    if vcf in "Varscan":
                        options.vcfformat = "vcf4old"
                    if vcf in "Lofreq":
                        options.exononly = False
                    success = 0
                    finalfn = vcf_filter(options)
                    finalout_dir = options.out_dir + "/finalResult"
                    create_dir(finalout_dir)
                    final_exon_fn = FundementalFile((str(finalfn)).replace("txt","exon.txt"))
                    final_exon_fn.cp("%s/%s.exon.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                    finalfn.cp("%s/%s.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                    success = 1
        if success:
            collect_result_file(vcfout_dir, options, bamfile_list, "germline", 2, 0.04)

def main():
    panel()


if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

