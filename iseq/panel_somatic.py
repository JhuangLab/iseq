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
    usage = "Usage: \n       %prog -c config.cfg -s A01 -m fastq2vcf -1 A01A_1.fq.gz -2 A01A_2.fq.gz -3 A01C_1.fq.gz -4 A01C_2.fq.gz --bamprocess 00101111 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -s A01 -m fastq2bam -1 A01A_1.fq -2 A01_2.fq -3 A01C_1.fq.gz -4 A01C_2.fq.gz --bamprocess 00101111 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -s A01 -m bam2vcf --case_in_bam A01A.bam --control_in_bam A01C.bam --bamprocess 00000000 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -m genomeindex \n"
    description = "iseq is an integrated analysis pipeline for NGS panel sequencing data. This is the control mode, case and control paired fastq and bam can be inputed. If you have any question about this tool, please contact us (lee_jianfeng@sjtu.edu.cn)"
    optparser = OptionParser(version = "0.1.0", description = description, usage = usage, add_help_option = False)
    optparser.add_option("-h", "--help", action = "help", help = "Show this help message and exit.")
    optparser.add_option("-c", "--config", dest = "config", default = "config.cfg" ,type = "string",
                         help = "Set the config File.[config.cfg]")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.(Required)")
    optparser.add_option("-m", "--mode", dest = "mode" ,type = "string",
                         help = "Run mode, [genomeindex, fastq2vcf, fastq2bam, bam2vcf, bamprocess, fastq2final, bam2final].(Required)")
    optparser.add_option("-1", "--case_fastq1", dest = "case_fastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-2", "--case_fastq2", dest = "case_fastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-3", "--control_fastq1", dest = "control_fastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-4", "--control_fastq2", dest = "control_fastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-d", "--bamprocess", dest = "bamprocess", default = "00101111",type = "string",
                         help = "Need point 8 digit number, eg. 00101111:  Add Read Group, Reorder Contig, Mark Duplicates, SplitNtrim ,RealignerTargetCreator,IndelRealigner,Recalibration,PrintReads, step one by one;01000000 only conduct Reorder contig; 00010000:Only conduct SplitNtrim step.[00101111].If --mode=Bamprocess,this parameter is required")
    optparser.add_option("-e", "--case_in_bam", dest = "case_in_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 and is empty, the in_bam is required!)")
    optparser.add_option("-b", "--control_in_bam", dest = "control_in_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 and is empty, the in_bam is required!)")
    optparser.add_option("-o", "--out_dir", dest = "out_dir" ,type = "string", default = "outdir",
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
    elif not isexist(options.config):
        optparser.print_help()
        print("Error:%s is not exist." % options.config)
        sys.exit(1)

    if options.mode is not None:
        options.mode = options.mode.lower()
    else:
        optparser.print_help()
        print("Error:Please set mode correctly.")
        sys.exit(1)

    if options.mode == "genomeindex":
        return(options)

    if not options.samplename:
        optparser.print_help()
        print("Error:Please set samplename correctly.")
        sys.exit(1)
    elif (not options.case_fastq1 or not options.control_fastq1) and (not options.case_in_bam or not options.control_in_bam):
        optparser.print_help()
        print("Error:Please set casefastq/controlfastq or case_in_bam/control_in_bam correctly.")
        sys.exit(1)
    elif options.mode not in ["genomeindex", "fastq2vcf", "fastq2bam", "bam2vcf", "bamprocess", "fastq2final", "bam2final"]:  
        optparser.print_help()
        print("Error:mode are not in genomeindex, fastq2vcf, fastq2bam, bam2vcf, bamprocess, fastq2final, bam2final.")
        sys.exit(1)
    elif options.mode in ["fastq2vcf","fastq2bam", "fastq2final"] and not options.case_fastq1 and not options.case_fastq2: 
        optparser.print_help()
        print("Error:Please set fastq1/fastq2 correctly." % options.config)
        sys.exit(1)
    elif options.mode in ["bam2vcf","bamprocess", "bam2final"] and not options.case_in_bam and not options.control_in_bam: 
        optparser.print_help()
        print("Error:Please set case_in_bam and control_in_bam correctly.")
        sys.exit(1)
    return(options)

def panel_somatic(options = ""):
    if options == "":
        options = opt_validate(prepare_optparser())
    if not isexist("restart"):
        os.makedirs("restart")
    create_dir("%s/log" % os.getcwd())
    cfg = get_config(options.config)
    vcf_out_dir = options.out_dir
    options.mode = options.mode.lower()
    ################ genome Index Mode #################
    if options.mode == "genomeindex":
        options.genome_index = "1"
        options.fastq_mapping = "0"
        options.bamprocess = "00000000"
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"

    ################ fastq2vcf Mode ################### 
    if options.mode == "fastq2vcf":
        options.genome_index = "0"
        options.fastq_mapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1"
        options.vcfannovar = "0"
        options.mpileup = "0"
    ################ fastq2final Mode ################### 
    if options.mode == "fastq2final":
        options.genome_index = "0"
        options.fastq_mapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1" 
        options.vcfannovar = "1"
        options.mpileup = "1"
    ################ fastq2bam Mode ################### 
    if options.mode == "fastq2bam":
        options.genomeindex = "0"
        options.fastq_mapping = "1"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"
        options.mpileup = "0"
    ################ bam2vcf Mode ################### 
    if options.mode == "bam2vcf":
        options.genome_index = "0"
        options.fastq_mapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1" 
        options.vcfannovar = "0"
        options.mpileup = "0"
    ################ bam2final Mode ################### 
    if options.mode == "bam2final":
        options.genome_index = "0"
        options.fastq_mapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "1"
        options.vcffilter = "1" 
        options.vcfannovar = "1"
        options.mpileup = "1"
    if options.mode == "bamprocess":
        options.genome_index = "0"
        options.fastq_mapping = "0"
        options.bamprocess = options.bamprocess
        options.variantcaller = "0"
        options.vcffilter = "0" 
        options.vcfannovar = "0"
        options.mpileup = "0"

    ################ Main Region ###################### 

    mapper = cfg["mapper"]
    mapper = mapper.split(",")
    frq_exon_only = cfg["freq_exon_only"]
    if frq_exon_only is "1":
      options.exononly = True
    else:
      options.exononly = False
    threads_mapper = []
    bamfiles_pool = {}
    if options.mode == "genomeindex":
       status = pre_process_somatic(options)
       return(status)
    else:
        for i in mapper:
            options.mapper = i
            def func(options = options):
                options = copy.deepcopy(options)
                case_bamfile_list,control_bamfile_list = pre_process_somatic(options)
                bamfiles_pool.update(case_bamfile_list)
                options.seq_type = "dna"
                if options.variantcaller == "1":
                    threads = []
                    for key in case_bamfile_list.keys():
                        def single(case_bamfile_list = case_bamfile_list, control_bamfile_list = control_bamfile_list, options = options, key = key):
                            options = copy.deepcopy(options)
                            case_bamfile = case_bamfile_list[key]
                            control_bamfile = control_bamfile_list[key]
                            options.in_bam = str(case_bamfile) + "," + str(control_bamfile)
                            options.out_dir = vcf_out_dir + "/call/" + key
                            options.runid = "%s.%s" % (options.samplename, key)
                            vcffile_list = variant_caller_somatic(options)
                            threads_vcf = []
                            for vcf in vcffile_list.keys():
                                def func_vcf(vcf = vcf, key = key, vcffile_list = vcffile_list, options = options):
                                    options = copy.deepcopy(options)
                                    key = key.capitalize()
                                    vcffile = vcffile_list[vcf]
                                    if vcffile:
                                        options.out_dir = vcffile.dirname + "/" 
                                        options.case_vcf = str(vcffile)
                                        vcf = vcf.capitalize()
                                    if vcf in "Unifiedgenotyper":
                                        options.vcffilter = cfg["unifiedgenotyper_filtration"]
                                    if vcf in "Varscan":
                                        options.vcfformat = "vcf4old"
                                        options.vcffilter = cfg["varscan_filtration"]
                                    if vcf in "Lofreq":
                                        options.vcffilter = cfg["lofreq_filtration"]
                                    if vcf in "Mutect":
                                        options.vcffilter = cfg["mutect_filtration"]
                                    if vcf in "Haplotypecaller":
                                        options.vcffilter = cfg["haplotypecaller_filtration"]
                                    options.runid = "%s.%s.%s" % (options.samplename, key, vcf)
                                    finalfn = vcf_filter_somatic(options)
                                    finalout_dir = options.out_dir + "/finalResult"
                                    if options.mode.find("final") != -1:
                                        create_dir(finalout_dir)
                                        final_exon_fn = FundementalFile((str(finalfn)).replace("txt","exon.txt"))
                                        final_exon_fn.cp("%s/%s.exon.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                                        finalfn.cp("%s/%s.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                                    else:
                                        finalfn.cp("%s/%s.%s.%s.vcf" % (finalout_dir, options.samplename, key, vcf))

                                threads_vcf.append(threading.Thread(target = func_vcf))
                            for t in threads_vcf:
                                t.setDaemon(True)
                                t.start()
                            for t in threads_vcf:
                                t.join()
                        threads.append(threading.Thread(target = single))
                    for t in threads:
                        t.setDaemon(True)
                        t.start()
                    for t in threads:
                        t.join()
            threads_mapper.append(threading.Thread(target = func))
        for t in threads_mapper:
            t.setDaemon(True)
            t.start()

        for t in threads_mapper:
            t.join()
        if options.mode.find("final") != -1:
            collect_result_file(vcf_out_dir, options, bamfiles_pool, "somatic",10, 0.068, 0.06)
        

def main():
    panel_somatic()


if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

