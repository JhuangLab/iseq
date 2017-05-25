#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2015 Jinyan HUANG <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License.
@Whole exon sequcing analysis pipeline in somatic mode
@status:  experimental
@version: 1.0 
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from DataClass import *
from PreProcess import *
from VariantCaller import *
from Refinement import *

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
                         help = "Run mode, [Genomeindex, Fastq2vcf, Fastq2bam, Bam2vcf, Bamprocess].(Required)")
    optparser.add_option("-1", "--casefastq1", dest = "casefastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-2", "--casefastq2", dest = "casefastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-3", "--controlfastq1", dest = "controlfastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-4", "--controlfastq2", dest = "controlfastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-d", "--bamprocess", dest = "bamprocess", default = "00101111",type = "string",
                         help = "Need point 8 digit number, eg. 00101111:  Add Read Group, Reorder Contig, Mark Duplicates, SplitNtrim ,RealignerTargetCreator,IndelRealigner,Recalibration,PrintReads, step one by one;01000000 only conduct Reorder contig; 00010000:Only conduct SplitNtrim step.[00101111].If --mode=Bamprocess,this parameter is required")
    optparser.add_option("-e", "--casein_bam", dest = "casein_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 and is empty, the in_bam is required!)")
    optparser.add_option("-b", "--controlin_bam", dest = "controlin_bam" ,type = "string",
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
    elif not options.casefastq1 and not options.controlfastq1 and not options.casein_bam and not options.controlin_bam:
        optparser.print_help()
        sys.exit(1)
    elif not options.mode:  
        optparser.print_help()
        sys.exit(1)
    elif options.mode not in ["Genomeindex", "Fastq2vcf", "Fastq2bam", "Bam2vcf", "Bamprocess"]:  
        optparser.print_help()
        sys.exit(1)
    elif options.mode in ["Fastq2vcf","Fastq2bam"] and not options.casefastq1 and not options.casefastq2: 
        optparser.print_help()
        sys.exit(1)
    elif options.mode in ["Bam2vcf","Bamprocess"] and not options.casein_bam and not options.controlin_bam: 
        optparser.print_help()
        sys.exit(1)
    return(options)

def weseqSomatic(options = ""):
    if options == "":
        options = opt_validate(prepare_optparser())
    cfg = get_config(options.config)
    vcf_out_dir = options.out_dir
    ################ Genome Index Mode #################
    if options.mode == "Genomeindex":
        options.Genomeindex = "1"
        options.FastqMapping = "0"
        options.bamprocess = "00000000"
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"

    ################ Fastq2vcf Mode ################### 
    if options.mode == "Fastq2vcf":
        options.Genomeindex = "0"
        options.FastqMapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1" 
        options.vcfannovar = "1"
        options.mpileup = "1"
    ################ Fastq2bam Mode ################### 
    if options.mode == "Fastq2bam":
        options.Genomeindex = "0"
        options.FastqMapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"
        options.mpileup = "0"
    ################ Bam2vcf Mode ################### 
    if options.mode == "Bam2vcf":
        options.Genomeindex = "0"
        options.FastqMapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1" 
        options.vcfannovar = "1"
        options.mpileup = "1"
    if options.mode == "Bamprocess":
        options.Genomeindex = "0"
        options.FastqMapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"
        options.mpileup = "0"

    ################ Main Region ###################### 
    casebamfileList,controlbamfileList = PreProcess_Somatic(options)
    options.seqtype = "dna"
    if options.variantcaller == "1":
        for key in casebamfileList.keys():
            casebamfile = casebamfileList[key]
            controlbamfile = controlbamfileList[key]
            options.in_bam = str(casebamfile) + "," + str(controlbamfile)
            options.out_dir = vcf_out_dir + "/" + key
            vcffileList = VariantCaller_Somatic(options)
            for vcf in vcffileList.keys():
                vcffile = vcffileList[vcf]
                if vcffile:
                    options.out_dir = vcffile.dirname + "/" 
                    options.case_vcf = str(vcffile)
                    if vcf in "UnifiedGenotyper":
                        pass
                        #options.vcffilter = "wespipeline"
                    options.exononly = True 
                    if vcf in "varscan":
                        options.vcfformat = "vcf4old"
                    if vcf in "lofreq":
                        options.exononly = False
                    success = 0
                    try:
                        finalfn = vcfFilter_somatic(options)
                        finalout_dir = options.out_dir + "/finalResult"
                        create_dir(finalout_dir)
                        finalexonfn = fundementalfile((str(finalfn)).replace("txt","exon.txt"))
                        finalexonfn.cp("%s/%s.exon.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                        finalfn.cp("%s/%s.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                        success = 1
                    except:
                        info("The " + vcf + " call result is empty!")
        if success:
            collect_result_file(vcf_out_dir, options.samplename, "somatic",10, 0.068, 0.06)
        

def main():
    weseqSomatic()


if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

