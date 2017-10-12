#! /usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@preprocess NGS Pre-Process Module
@status:  experimental
@version: 0.0.1
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from fastq import FastqFile
from reffa import ReffaFile
from bam import BamFile

def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = "usage: %prog -c mysample.cfg -s A01A -1 A01_1.fq -2 A02_2.fq"
    description = "Please set the sample name. e.g. L04A, L04C, L04T."
    optparser = OptionParser(version = "0.0.1", description = description, usage = usage, add_help_option = False)
    optparser.add_option("-h", "--help", action = "help", help = "Show this help message and exit.")
    optparser.add_option("-c", "--config", dest = "config", default = "config.cfg" ,type = "string",
                         help = "Set the config File.")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.")
    optparser.add_option("-1", "--fastq1", dest = "fastq1", type = "string",
                         help = "input fastq file paired 1")
    optparser.add_option("-2", "--fastq2", dest = "fastq2", type = "string",
                         help = "input fastq file paired 2")
    optparser.add_option("-d", "--dataprocess", dest = "dataprocess", default = "1111111111",type = "string",
                         help = "Need point 6 digit number, eg. 111111: Conduct Genome Process, fastq_mapping, Add Read Group, Reorder Contig, Mark Duplicates, split_ntrim step one by one;100000 only conduct Genome Process; 000001:Only conduct split_ntrim step")
    optparser.add_option("-i", "--in_bam", dest = "in_bam" ,type = "string",
                         help = "You can set this to your bam file path.(If fastq1 is empty, required!)")
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
    elif not options.fastq1 and not options.fastq2 and not options.in_bam:
        optparser.print_help()
        sys.exit(1)
    elif options.dataprocess[2] == "1" and not options.fastq1 and not options.fastq2:
        optparser.print_help()
        sys.exit(1)
    return(options)
class FundementalPreprocess(object):
    def __init__(self, options):
        self.samplename = options.samplename
        self.options = options
        self.cfg = get_config(options.config)
        self.reffa = self.cfg["reffa"]
        try:
            self.mapper = options.mapper
        except:
            pass
        self.genome_indexer = self.cfg["genome_indexer"]

class GenomePreprocessor(FundementalPreprocess):
    def __init__(self, options):
        FundementalPreprocess.__init__(self, options)
        self.genome_indexer_list = self.genome_indexer.split(",") 
    def preprocess(self):
        genome = ReffaFile(self.reffa, self.cfg)
        threads = []
        for indexer in self.genome_indexer_list:
            indexer = indexer.lower()
            if indexer == "star":
                #1pass to get SJ.tab file for run 2pass when run STAR mapping step
                t1 = threading.Thread(target = genome.star_index)
                threads.append(t1)
            elif indexer == "bwa":
                t2 = threading.Thread(target = genome.bwa_index)
                threads.append(t2)
            elif indexer == "bowtie2":
                t3 = threading.Thread(target = genome.bowtie2_index)
                threads.append(t3)
            elif indexer == "bowtie":
                t4 = threading.Thread(target = genome.bowtie_index)
                threads.append(t4)
            elif indexer == "tmap":
                t4 = threading.Thread(target = genome.tmap_index)
                threads.append(t4)
            else:
                info("The genome_indexer set error in config.cfg!")

        t5 = threading.Thread(target = genome.generate_dict)
        threads.append(t5)
        # ReffaFile Class instance
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(genome)

class FastqPreprocessor(FundementalPreprocess):
    def __init__(self, options):
        FundementalPreprocess.__init__(self, options)
        self.starbam = self.options.out_dir + "/bam/Star" 
        self.bwabam = self.options.out_dir + "/bam/Bwa"
        self.bowtie2bam = self.options.out_dir + "/bam/Bowtie2"
        self.bowtiebam = self.options.out_dir + "/bam/Bowtie"
        self.tophatbam = self.options.out_dir + "/bam/Tophat"
        self.tophatbam = self.options.out_dir + "/bam/Tmap"
        self.fastq1 = options.fastq1
        self.fastq2 = options.fastq2
        self.fastq1fn = FastqFile(self.fastq1 ,self.samplename, self.cfg)
        self.fastq2fn = self.fastq2
        #Use config mapper parameter to conduct mapping step one by one
        self.mapper_list = self.mapper.split(",") 
        self.bamfile_list = {}
    def preprocess(self):
        """
        Ret:Need two or one fastq file, and use config file point mapper to conduct mapping step.After run this method, it will set the self.fastq_mappingBamList value.
        """
        genome = ReffaFile(self.reffa, self.cfg)
        threads = []
        for mapper in self.mapper_list:
            mapper = mapper.upper()
            if mapper == "STAR":  # STAR mapping step need two pass, 2pass use 1pass generating SJ.out.tab to generate Genome index
                def star_single (self = self):
                    out_bam_dir =  "%s/%s/1pass/" % (self.starbam, self.samplename)
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.star_mapping(out_bam_dir, self.fastq2fn)
                    out_bam_dir =  "%s/%s/2pass/" % (self.starbam, self.samplename)
                    create_dir(out_bam_dir)
                    genome_2pass_out_dir = "%s/%s/genome_2pass" % (self.starbam, self.samplename)
                    genome.runid = self.fastq1fn.runid + ".genome_2pass"
                    genome.star_index("2pass", genome_2pass_out_dir, self.samplename, out_bam.dirname + "/SJ.out.tab")
                    reffa = genome_2pass_out_dir + "/genome" 
                    cfg = copy.deepcopy(self.cfg)
                    cfg["reffa"] = reffa
                    self.fastq1fn.runid = self.fastq1fn.runid + ".2pass"
                    out_bam = self.fastq1fn.star_mapping(out_bam_dir, self.fastq2fn)
                    #If out_bam is not be generated, it can be bool value False, So need to be try/excpet to eat warning msg
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Star"] = out_bam 
                    except:
                        info(str(self.fastq1fn) + " STAR mapping step run unfinished!")
                t1 = threading.Thread(target = star_single)
                threads.append(t1)
            elif mapper == "BWA":
                def bwa_single(self=self):
                    out_bam_dir =  self.bwabam + "/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.bwa_mapping(out_bam_dir, self.fastq2fn)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Bwa"] = out_bam 
                    except:
                        info(str(self.fastq1fn) + " Bwa mapping step run unfinished!")
                t2 = threading.Thread(target = bwa_single)
                threads.append(t2)
            elif mapper == "BOWTIE2":
                def bowtie2_single(self=self):
                    out_bam_dir =  self.bowtie2bam + "/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.bowtie2_mapping(out_bam_dir, self.fastq2fn)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Bowtie2"] = out_bam 
                    except:
                        info(str(self.fastq1fn) + " Bwa mapping step run unfinished!")
                t3 = threading.Thread(target = bowtie2_single)
                threads.append(t3)
            elif mapper == "BOWTIE":
                def bowtie_single(self=self):
                    out_bam_dir =  self.bowtiebam + "/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.bowtie_mapping(out_bam_dir, self.fastq2fn)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Bowtie"] = out_bam 
                    except:
                        info(str(self.fastq1fn) + " Bowtie mapping step run unfinished!")
                t4 = threading.Thread(target = bowtie_single)
                threads.append(t4)
            elif mapper == "TOPHAT2":
                def tophat2_single(self=self):
                    out_bam_dir =  self.tophatbam + "2/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.tophat_mapping(out_bam_dir, "tophat2", self.fastq2fn)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Tophat2"] = out_bam
                    except:
                        info(str(self.fastq1fn) + " Tophat mapping using bowtie2 step run unfinished!")
                t5 = threading.Thread(target = tophat2_single)
                threads.append(t5)
            elif mapper == "TOPHAT":
                def tophat_single(self=self):
                    out_bam_dir =  self.tophatbam + "/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.tophat_mapping(out_bam_dir, "tophat", self.fastq2fn)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Tophat"] = out_bam
                    except:
                        info(str(self.fastq1fn) + " Tophat mapping using bowtie step run unfinished!")
                t6 = threading.Thread(target = tophat2_single)
                threads.append(t6)
            elif mapper == "TMAP":
                def tmap_single(self=self):
                    out_bam_dir =  self.tmapbam + "/" + self.samplename
                    create_dir(out_bam_dir)
                    out_bam = self.fastq1fn.tmap_mapping(out_bam_dir)
                    try:
                        if out_bam.isexist():
                            self.bamfile_list["Tmap"] = out_bam
                    except:
                        info(str(self.fastq1fn) + " TMAP mapping step run unfinished!")
                t7 = threading.Thread(target = tmap_single)
                threads.append(t7)
            else:
                info("The mapper set error in config.cfg!")
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list) # the dict of BamFile class instance 


class BamPreprocessor(FundementalPreprocess):
    """
    Ret: The class is be used to process bam file to final can be conduct caller bam
    """
    def __init__(self, options, bamfile_list):
        FundementalPreprocess.__init__(self, options)
        self.bamfile_list = bamfile_list
    def add_read_group(self):
        """
        Ret:Need set the self.bamfile_list value, it will use picard to set  RGID = 1, RGLB = "Jhuanglab", RGPL="ILLUMINA", RGPU = "Hiseq",RGSM=self.samplename
        """
        threads = []
        for key in self.bamfile_list.keys():
            def add_read_group_single(self = self):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_AddGroup.bam"
                out_bam, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                groupedBam = bamfile.add_read_group(out_bam)
                #groupedBam.index(self.samtools)
                self.bamfile_list[key] = groupedBam
            threads.append(threading.Thread(target = add_read_group_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def contig_reorder(self):
        """
        Ret:Need set the self.bamfile_list value, it will use picard to reoder bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def contig_reorder_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_Reorder.bam"
                out_bam, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                reorderedBam = bamfile.contig_reorder(out_bam)
                #reorderedBam.index(self.samtools)
                self.bamfile_list[key] = reorderedBam 
            threads.append(threading.Thread(target = contig_reorder_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def mark_duplicates(self):
        """
        Ret:Need set the self.bamfile_list value, it will use picard to mark duplcates with bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def mark_duplicates_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_MarkDup.bam"
                out_bam, number = re.subn(pattern, new_string, bamfile.path)
                new_string = ".metrics"
                out_metrics, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                mark_duplicatesBam = bamfile.mark_duplicates(out_bam, out_metrics)
                #mark_duplicatesBam.index(self.samtools)
                self.bamfile_list[key] = mark_duplicatesBam 
            threads.append(threading.Thread(target = mark_duplicates_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def realigner_target_creator(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK realigner_target_creator to bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def realigner_target_creator_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = ".intervals"
                out_intervals, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                bamfile.realigner_target_creator(out_intervals)
                bamfile.intervals = out_intervals
                self.bamfile_list[key] = bamfile 
            threads.append(threading.Thread(target = realigner_target_creator_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def indel_realigner(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK indel_realigner to bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def indel_realigner_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_IndelRealigner.bam"
                out_bam, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                indel_realignerBam = bamfile.indel_realigner(
                        bamfile.intervals, out_bam)
                self.bamfile_list[key] = indel_realignerBam 
            threads.append(threading.Thread(target = indel_realigner_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def recalibration(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK BaseRecalibrator to bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def recalibration_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_Recal_data.grp"
                out_grp, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                bamfile.recalibration(out_grp)
                bamfile.recaldata = out_grp
                self.bamfile_list[key] = bamfile 
            threads.append(threading.Thread(target = recalibration_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def print_reads(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK to generate preprocessed bam File in self.bamfile_list
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def print_reads_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = ".bam$"
                new_string = "_PrintReads.bam"
                out_bam, number = re.subn(pattern, new_string, bamfile.path)
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                print_readsBam = bamfile.print_reads(bamfile.recaldata, out_bam)
                self.bamfile_list[key] = print_readsBam 
            threads.append(threading.Thread(target = print_reads_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)
    def split_ntrim(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK to split_ntrim with bam File in self.bamfile_list, it will set the self.split_ntrimBam value.
        """ 
        threads = []
        for key in self.bamfile_list.keys():
            def split_ntrim_single(self = self, key = key):
                bamfile = self.bamfile_list[key]
                pattern = re.compile(".bam$")
                obj = re.search(pattern, bamfile.path)
                replace_str = obj.group()
                out_bam = bamfile.path.replace(replace_str, "_SplitNtrim.bam")
                if key not in bamfile.runid:
                    bamfile.runid = bamfile.runid + "." + key
                splitNtrimBam = bamfile.split_ntrim(out_bam)
                splitNtrimBam.index()
                self.bamfile_list[key] = splitNtrimBam
            threads.append(threading.Thread(target = split_ntrim_single))
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.bamfile_list)



class PreProcessor(FundementalPreprocess):
    """
    Description:
        Analysis Pre-Process Class,need a run option object initial. 
    Method:
        genome_pre_process:Use STAR/bwa/bowtie2/bowtie to index Reference.
        fastq_mapping:Need two or one fastq file, and use config file point mapper to conduct mapping step.After run this method, it will set the self.fastq_mappingBamList value. 
        add_read_group:Need set the self.bamfile_list value, it will use picard to set  RGID = 1, RGLB = "Jhuanglab", RGPL="ILLUMINA", RGPU = "Hiseq",RGSM=self.samplename, it will set the self.add_read_groupBamList value. 
        contig_reorder:Need set the self.bamfile_list value, it will use picard to reoder bam File in self.bamfile_list, it will set the self.contig_reorderBamList value.
        mark_duplicates:Need set the self.bamfile_list value, it will use picard to mark duplcates with bam File in self.bamfile_list, it will set the self.mark_duplicatesBamList value.
        split_ntrim:Need set the self.bamfile_list value, it will use GATK to split_ntrim with bam File in self.bamfile_list, it will set the self.split_ntrimBam value.
    """
    def __init__(self, options):
        FundementalPreprocess.__init__(self, options)
        if options.mode != "genomeindex":
            self.bamfile_list = {"Default":BamFile(options.in_bam, self.samplename, self.cfg) } 

    def genome_pre_process(self):
        """
        Ret:Use STAR/bwa/bowtie2/bowtie to index Reference.
        """
        Genomeprocess = GenomePreprocessor(self.options) 
        Genomeprocess.preprocess()

    def fastq_mapping(self):
        Fastqmapping = FastqPreprocessor(self.options) 
        self.bamfile_list = Fastqmapping.preprocess()
        self.fastq_mappingBamList = self.bamfile_list 
        return(self.bamfile_list)
    def add_read_group(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.add_read_group()
        self.add_read_groupBamList = self.bamfile_list 
        return(self.bamfile_list)
    def contig_reorder(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.contig_reorder()
        self.contig_reorderBamList = self.bamfile_list 
        return(self.bamfile_list)
    def mark_duplicates(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.mark_duplicates()
        self.mark_duplicatesBamList = self.bamfile_list 
        return(self.bamfile_list)
    def realigner_target_creator(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.realigner_target_creator()
        self.realigner_target_creatorBamList = self.bamfile_list 
        return(self.bamfile_list)
    def indel_realigner(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.indel_realigner()
        self.indel_realignerBamList = self.bamfile_list 
        return(self.bamfile_list)
    def recalibration(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.recalibration()
        self.recalibrationBamList = self.bamfile_list 
        return(self.bamfile_list)
    def print_reads(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.print_reads()
        self.print_readsBamList = self.bamfile_list 
        return(self.bamfile_list)
    def split_ntrim(self):
        bamfiles = BamPreprocessor(self.options, self.bamfile_list)
        self.bamfile_list = bamfiles.split_ntrim()
        self.split_ntrimBam = self.bamfile_list 
        return(self.bamfile_list)
class PreProcessorSomatic(FundementalPreprocess):
    """
    Description:
        Analysis Pre-Process Class,need a run option object initial. 
    Method:
        genome_pre_process:Use STAR/bwa/bowtie2/bowtie to index Reference.
        fastq_mapping:Need two or one fastq file, and use config file point mapper to conduct mapping step.After run this method, it will set the self.fastq_mappingBamList value. 
        add_read_group:Need set the self.bamfile_list value, it will use picard to set  RGID = 1, RGLB = "Jhuanglab", RGPL="ILLUMINA", RGPU = "Hiseq",RGSM=self.samplename, it will set the self.add_read_groupBamList value. 
        contig_reorder:Need set the self.bamfile_list value, it will use picard to reoder bam File in self.bamfile_list, it will set the self.contig_reorderBamList value.
        mark_duplicates:Need set the self.bamfile_list value, it will use picard to mark duplcates with bam File in self.bamfile_list, it will set the self.mark_duplicatesBamList value.
        split_ntrim:Need set the self.bamfile_list value, it will use GATK to split_ntrim with bam File in self.bamfile_list, it will set the self.split_ntrimBam value.
    """
    def __init__(self, options):
        FundementalPreprocess.__init__(self, options)
        if options.mode != "genomeindex":
            self.case_bamfile_list = {"Bwa":BamFile(options.case_in_bam, self.samplename + "A", self.cfg)}
            self.control_bamfile_list = {"Bwa":BamFile(options.control_in_bam, self.samplename + "C", self.cfg)}
            self.caseoptions = copy.deepcopy(options)
            self.caseoptions.fastq1 = self.options.case_fastq1 
            self.caseoptions.fastq2 = self.options.case_fastq2 
            self.caseoptions.samplename = options.samplename + "A"
            self.controloptions = copy.deepcopy(options)
            self.controloptions.fastq1 = options.control_fastq1 
            self.controloptions.fastq2 = options.control_fastq2 
            self.controloptions.samplename = options.samplename + "C"

    def genome_pre_process(self):
        """
        Ret:Use STAR/bwa/bowtie2/bowtie to index Reference.
        """
        Genomeprocess = GenomePreprocessor(self.options) 
        Genomeprocess.preprocess()

    def fastq_mapping(self):
        def case_func(self = self):
            #case
            caseoptions = self.caseoptions
            caseFastqmapping = FastqPreprocessor(caseoptions) 
            self.case_bamfile_list = caseFastqmapping.preprocess()
            self.casefastq_mappingBamList = self.case_bamfile_list 
        def control_func(self = self):
            #control
            controloptions = self.controloptions
            controlFastqmapping = FastqPreprocessor(controloptions) 
            self.control_bamfile_list = controlFastqmapping.preprocess()
            self.controlfastq_mappingBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def add_read_group(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.add_read_group()
            self.caseadd_read_groupBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.add_read_group()
            self.controladd_read_groupBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def contig_reorder(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.contig_reorder()
            self.casecontig_reorderBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.contig_reorder()
            self.controlcontig_reorderBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def mark_duplicates(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.mark_duplicates()
            self.casemark_duplicatesBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.mark_duplicates()
            self.controlmark_duplicatesBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def realigner_target_creator(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.realigner_target_creator()
            self.caserealigner_target_creatorBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.realigner_target_creator()
            self.controlrealigner_target_creatorBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def indel_realigner(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.indel_realigner()
            self.caseindel_realignerBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.indel_realigner()
            self.controlindel_realignerBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def recalibration(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.recalibration()
            self.caserecalibrationBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.recalibration()
            self.controlrecalibrationBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def print_reads(self):
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.print_reads()
            self.caseprint_readsBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.print_reads()
            self.controlprint_readsBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)
    def split_ntrim(self):
        """
        Ret:Need set the self.bamfile_list value, it will use GATK to split_ntrim with bam File in self.bamfile_list, it will set the self.split_ntrimBam value.
        """ 
        def case_func(self = self):
            case_bamfiles = BamPreprocessor(self.caseoptions, self.case_bamfile_list)
            self.case_bamfile_list = case_bamfiles.split_ntrim()
            self.casesplit_ntrimBamList = self.case_bamfile_list 
        def control_func(self = self):
            control_bamfiles = BamPreprocessor(self.controloptions, self.control_bamfile_list)
            self.control_bamfile_list = control_bamfiles.split_ntrim()
            self.controlsplit_ntrimBamList = self.control_bamfile_list 
        threads = []
        t1 = threading.Thread(target = case_func) 
        t2 = threading.Thread(target = control_func) 
        threads.append(t1)
        threads.append(t2)
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
        return(self.case_bamfile_list,self.control_bamfile_list)


def pre_process(options):
    """
    pre_process:Be imported in rnaseq.py, conduct the genome index, fastq mapping, several of bam file preprocess step.
    """
    cfg = get_config(options.config)
    preprocess =  PreProcessor(options)
    try:
        genome_process = options.genome_index #Genomeindex be set in rnaseq.py mode parameter
    except:
        genome_process = 0
    try:
        fastq_mapping = options.fastq_mapping #fastq_mapping be set in rnaseq.py mode parameter
    except:
        fastq_mapping = 0
    add_read_group = options.bamprocess[0]    # -d --bamprocess 0111
    contig_reorder = options.bamprocess[1]   # 0 or 1
    mark_duplicates = options.bamprocess[2]
    split_ntrim = options.bamprocess[3] 
    realigner_target_creator = options.bamprocess[4] 
    indel_realigner = options.bamprocess[5] 
    recalibration = options.bamprocess[6] 
    print_reads = options.bamprocess[7] 
    ############# Value equal 0 not to be run, and equal 1 will be run ##########
    if genome_process != "0":
        status = preprocess.genome_pre_process()
    if options.mode == "genomeindex":
        return(status)
    if fastq_mapping != "0":
        preprocess.fastq_mapping()
    #Bam process
    if add_read_group != "0":
        preprocess.add_read_group()
    if contig_reorder != "0":
        preprocess.contig_reorder()
    if mark_duplicates != "0":
        preprocess.mark_duplicates()
    if split_ntrim != "0":
        preprocess.split_ntrim()
    if realigner_target_creator != "0":
        preprocess.realigner_target_creator()
    if indel_realigner != "0":
        preprocess.indel_realigner()
    if recalibration != "0":
        preprocess.recalibration()
    if print_reads != "0":
        preprocess.print_reads()
    #Get copy and rename
    for key in preprocess.bamfile_list.keys():
        bamfile = preprocess.bamfile_list[key]
        processed_bam =  bamfile#.dirname + "/" + bamfile.samplename+"_preprocessed.bam"
        if bamfile.isexist(): #and bamfile.ln(processed_bam, bamfile.samplename):
            bamfile = BamFile(processed_bam, bamfile.samplename, cfg)
            #bamfile.index(cfg["samtools"])
            preprocess.bamfile_list[key] = bamfile
    return(preprocess.bamfile_list) 

def pre_process_somatic(options):
    """
    pre_process_somatic:Be imported in somatic.py, conduct the genome index, fastq mapping, several of bam file preprocess step.
    """
    cfg = get_config(options.config)
    preprocess =  PreProcessorSomatic(options)
    try:
        genome_process = options.genome_index #Genomeindex be set in somatic.py mode parameter
    except:
        genome_process = 0
    try:
        fastq_mapping = options.fastq_mapping #fastq_mapping be set in somatic.py mode parameter
    except:
        fastq_mapping = 0
    add_read_group = options.bamprocess[0]    # -d --bamprocess 0111
    contig_reorder = options.bamprocess[1]   # 0 or 1
    mark_duplicates = options.bamprocess[2]
    split_ntrim = options.bamprocess[3] 
    realigner_target_creator = options.bamprocess[4] 
    indel_realigner = options.bamprocess[5] 
    recalibration = options.bamprocess[6] 
    print_reads = options.bamprocess[7] 
    ############# Value equal 0 not to be run, and equal 1 will be run ##########
    if genome_process != "0":
        status = preprocess.genome_pre_process()
    if options.mode == "genomeindex":
        return(status)

    if fastq_mapping != "0":
        preprocess.fastq_mapping()
    #Bam process
    if add_read_group != "0":
        preprocess.add_read_group()
    if contig_reorder != "0":
        preprocess.contig_reorder()
    if mark_duplicates != "0":
        preprocess.mark_duplicates()
    if split_ntrim != "0":
        preprocess.split_ntrim()
    if realigner_target_creator != "0":
        preprocess.realigner_target_creator()
    if indel_realigner != "0":
        preprocess.indel_realigner()
    if recalibration != "0":
        preprocess.recalibration()
    if print_reads != "0":
        preprocess.print_reads()
    #Get copy and rename
    threads = []
    for key in preprocess.case_bamfile_list.keys():
        def case_func(preprocess = preprocess, key = key):
            case_bamfile = preprocess.case_bamfile_list[key]
            processed_bam =  case_bamfile#case_bamfile.dirname + "/" + case_bamfile.samplename#+"_preprocessed.bam"
            if case_bamfile.isexist(): #and case_bamfile.ln(processed_bam, case_bamfile.samplename):
                case_bamfile = BamFile(processed_bam, case_bamfile.samplename, cfg)
                #case_bamfile.index(cfg["samtools"])
                preprocess.case_bamfile_list[key] = case_bamfile
        def control_func(preprocess = preprocess, key =key):
            control_bamfile = preprocess.control_bamfile_list[key]
            processed_bam =  control_bamfile#control_bamfile.dirname + "/" + control_bamfile.samplename#+"_preprocessed.bam"
            if control_bamfile.isexist(): #and control_bamfile.ln(processed_bam, control_bamfile.samplename):
                control_bamfile = BamFile(processed_bam, control_bamfile.samplename, cfg)
                #control_bamfile.index(cfg["samtools"])
                preprocess.control_bamfile_list[key] = control_bamfile
        threads.append(threading.Thread(target = case_func))
        threads.append(threading.Thread(target = control_func)) 
        for t in threads:
            t.setDaemon(True)
            t.start()

        for t in threads:
            t.join()
    return(preprocess.case_bamfile_list, preprocess.control_bamfile_list)

def pre_process_self(options):
    """
    pre_process_self:Be useed, when run the module of self, using options.dataprocess to control the step be runned.
    """
    cfg = get_config(options.config)
    preprocess =  PreProcessor(options)
    genome_process = options.dataprocess[0]
    fastq_mapping = options.dataprocess[1]
    add_read_group = options.dataprocess[2]
    contig_reorder = options.dataprocess[3]
    mark_duplicates = options.dataprocess[4]
    split_ntrim = options.dataprocess[5] 
    realigner_target_creator = options.dataprocess[6] 
    indel_realigner = options.dataprocess[7] 
    recalibration = options.dataprocess[8] 
    print_reads = options.dataprocess[9] 
    if genome_process != "0":
        status = preprocess.genome_pre_process()
    if options.mode == "genomeindex":
        return(status)
    if fastq_mapping != "0":
        preprocess.fastq_mapping()
    if add_read_group != "0":
        preprocess.add_read_group()
    if contig_reorder != "0":
        preprocess.contig_reorder()
    if mark_duplicates != "0":
        preprocess.mark_duplicates()
    if split_ntrim != "0":
        preprocess.split_ntrim()
    if realigner_target_creator != "0":
        preprocess.realigner_target_creator()
    if indel_realigner != "0":
        preprocess.indel_realigner()
    if recalibration != "0":
        preprocess.recalibration()
    if print_reads != "0":
        preprocess.print_reads()
    for key in preprocess.bamfile_list.keys():
        bamfile = preprocess.bamfile_list[key]
        processed_bam =  bamfile#bamfile.dirname + "/" + bamfile.samplename#+"_preprocessed.bam"
        if bamfile.isexist(): #and bamfile.ln(processed_bam, bamfile.samplename):
            bamfile = BamFile(processed_bam, bamfile.samplename, cfg)
            #bamfile.index(cfg["samtools"])
            preprocess.bamfile_list[key] = bamfile
        else:
            info("Not found avaliable bam file in input!")
            return(False)

    return(preprocess.bamfile_list)

def main():
    options = opt_validate(prepare_optparser())
    bamfile_list = pre_process_self(options)
        

if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

