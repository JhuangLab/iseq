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
    usage = "Usage: \n       %prog -c config.cfg -s A01A -m fastq2vcf -1 A01_1.fq -2 A02_2.fq --bamprocess 00101111 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -s A01A -m fastq2bam -1 A01_1.fq -2 A02_2.fq --bamprocess 00101111 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -s A01A -m bam2vcf --in_bam A01.bam --bamprocess 00000000 -o outdir \n"
    usage = usage + "       %prog -c config.cfg -m genomeindex \n"
    description = "iseq is an integrated analysis pipeline for NGS panel sequencing data. This is the no control mode, fastq and bam can be inputed. If you have any question about this tool, please contact us (lee_jianfeng@sjtu.edu.cn)"
    optparser = OptionParser(version = "0.1.0", description = description, usage = usage, add_help_option = False)
    optparser.add_option("-h", "--help", action = "help", help = "Show this help message and exit.")
    optparser.add_option("-c", "--config", dest = "config", default = "config.cfg" ,type = "string",
                         help = "Set the config File.[config.cfg]")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.(Required)")
    optparser.add_option("-m", "--mode", dest = "mode" ,type = "string",
                         help = "Run mode, [genome_index, fastq2vcf, fastq2bam, bam2vcf, bamprocess, fastq2final, bam2final].(Required)")
    optparser.add_option("-1", "--fastq1", dest = "fastq1", type = "string", default = "",
                         help = "input fastq file paired 1.")
    optparser.add_option("-2", "--fastq2", dest = "fastq2", type = "string", default = "",
                         help = "input fastq file paired 2.")
    optparser.add_option("-d", "--bamprocess", dest = "bamprocess", default = "00101111",type = "string",
                         help = "Need point 8 digit number, eg. 00101111:  Add Read Group, Reorder Contig, Mark Duplicates, SplitNtrim ,RealignerTargetCreator,IndelRealigner,Recalibration,PrintReads, step one by one;01000000 only conduct Reorder contig; 00010000:Only conduct SplitNtrim step.[00101111].If --mode=Bamprocess,this parameter is required")
    optparser.add_option("-i", "--in_bam", dest = "in_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 and is empty, the in_bam is required!)")
    optparser.add_option("-o", "--out_dir", dest = "out_dir" ,type = "string", default = "outdir",
                         help = "Set the vcf file out_dir.[outdir]")
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
    elif not options.fastq1 and not options.fastq2 and not options.in_bam:
        optparser.print_help()
        print("Error:Please set fastq1/fastq2 or in_bam correctly.")
        sys.exit(1)
    elif options.mode not in ["genomeindex", "fastq2vcf", "fastq2bam", "bam2vcf", "bamprocess", "fastq2final", "bam2final"]:  
        optparser.print_help()
        print("Error:mode are not in genomeindex, fastq2vcf, fastq2bam, bam2vcf, bamprocess.")
        sys.exit(1)
    elif options.mode in ["fastq2vcf","fastq2bam", "fastq2final"] and not options.fastq1 and not options.fastq2: 
        optparser.print_help()
        print("Error:Please set fastq1/fastq2 correctly.")
        sys.exit(1)
    elif options.mode in ["bam2vcf","bamprocess", "bam2final"] and not options.in_bam: 
        optparser.print_help()
        print("Error:Please set in_bam correctly.")
        sys.exit(1)
    return(options)

def panel(options=""):
    if options == "":
        options = opt_validate(prepare_optparser())
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
        options.genome_index = "0"
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
    ################ bam2vcf Mode ################### 
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
        options.vcffilter = "1" 
        options.vcfannovar = "0"
        options.mpileup = "0"

    ################ Main Region ###################### 
    mapper = cfg["mapper"]
    mapper = mapper.split(",")
    threads_mapper = []
    bamfiles_pool = {}
    if options.mode == "genomeindex":
       status = pre_process(options)
       return(status)
    else:
        for i in mapper:
            options = copy.deepcopy(options)
            options.mapper = i
            def func(options = options, bamfiles_pool = bamfiles_pool):
                options = copy.deepcopy(options)
                bamfile_list = pre_process(options)
                bamfiles_pool.update(bamfile_list)
                options.seq_type = "dna"
                if options.variantcaller == "1":
                    threads_bam = []
                    for key in bamfile_list.keys():
                        def single(bamfile_list = bamfile_list, options = options, key = key):
                            options = copy.deepcopy(options)
                            bamfile = bamfile_list[key]
                            options.in_bam = str(bamfile)
                            options.out_dir = vcf_out_dir + "/call/" + key
                            options.runid = "%s.%s" % (options.samplename, key)
                            vcffile_list = variant_caller(options)
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
                                    options.exononly = True 
                                    if vcf in "Unifiedgenotyper":
                                        options.vcffilter = cfg["unifiedgenotyper_filtration"]
                                    if vcf in "Varscan":
                                        options.vcfformat = "vcf4old"
                                        options.vcffilter = cfg["varscan_filtration"]
                                    if vcf in "Lofreq":
                                        options.exononly = False
                                        options.vcffilter = cfg["lofreq_filtration"]
                                    if vcf in "Mutect":
                                        options.vcffilter = cfg["mutect_filtration"]
                                    if vcf in "Haplotypecaller":
                                        options.vcffilter = cfg["haplotypecaller_filtration"]
                                    options.runid = "%s.%s.%s" % (options.samplename, key, vcf)
                                    finalfn = vcf_filter(options)
                                    finalout_dir = options.out_dir + "/finalResult"
                                    create_dir(finalout_dir)
                                    final_exon_fn = FundementalFile((str(finalfn)).replace("txt","exon.txt"))
                                    final_exon_fn.cp("%s/%s.exon.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                                    finalfn.cp("%s/%s.%s.%s.txt" % (finalout_dir, options.samplename, key, vcf))
                                threads_vcf.append(threading.Thread(target = func_vcf))
                            for t in threads_vcf:
                                t.setDaemon(True)
                                t.start()
                            for t in threads_vcf:
                                t.join()
                        threads_bam.append(threading.Thread(target = single))
                    for t in threads_bam:
                        t.setDaemon(True)
                        t.start()
                    for t in threads_bam:
                        t.join()
            threads_mapper.append(threading.Thread(target = func))
        for t in threads_mapper:
            t.setDaemon(True)
            t.start()
        count = 0
        for t in threads_mapper:
            t.join()
            info("%s %s full process be finished, wating another mapper or continue.", options.samplename, mapper[count])
            count = count + 1
        if options.mode.find("final") != -1:
            collect_result_file(vcf_out_dir, options, bamfiles_pool, "germline", 2, 0.04)

def main():
    create_dir("%s/log" % os.getcwd())
    panel()


if __name__ == "__main__":
    try:
        create_dir("%s/log" % os.getcwd())
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

