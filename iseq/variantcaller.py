#! /usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li<lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@variantcaller Variant Call Module
@status: experimental
@version: $Revision$
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from bam import *
from vcf import *

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
    usage = "usage: %prog -c mysample.cfg -s A01A -1 A01_1.fq -2 A02_2.fq"
    description = "Please set the sample name. e.g. L04A, L04C, L04T."
    optparser = OptionParser(version = "0.0.1", description = description, usage = usage, add_help_option = False)
    optparser.add_option("-h", "--help", action = "help", help = "Show this help message and exit.")
    optparser.add_option("-c", "--config", dest = "config", default = "config.cfg" ,type = "string",
                         help = "Set the config File.")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.")
    optparser.add_option("-i", "--in_bam", dest = "in_bam" ,type = "string",
                         help = "If you have process the preprocess step, you can set this to your bam file path.(eg. case,control, or case only)")
    optparser.add_option("-t", "--seq_type", dest = "seq_type" ,type = "string", default = "dna",
                         help = "Point the seq type[dna]")
    optparser.add_option("-o", "--out_dir", dest = "out_dir" ,type = "string",
                         help = "Set the vcf file out_dir.(Required)")
    return(optparser)
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options, args) = optparser.parse_args()
    if not options.config:
        optparser.print_help()
        sys.exit(1)
    elif not options.in_bam:
        optparser.print_help()
        sys.exit(1)
    elif not options.out_dir:
        optparser.print_help()
        sys.exit(1)
    return(options)
class FundementalCaller(object):
    def __init__(self, caller_name):
        self.mode = "germline"
        self.avaliable_caller = ["lofreq","varscan","mutect","tvc","haplotypecaller","unifiedgenotyper","pindel"]
        if caller_name.lower() in self.avaliable_caller:
            self.caller_name = caller_name.lower()
        else:
            info("Now the avaliable variant_caller is only:" + self.avaliable_caller.join(",") + "! Please set correct caller in config file.")
    def set_out_dir(self, out_dir):
        self.out_dir = out_dir 
        create_dir(out_dir)
    def set_caller_mode(self, mode):
        if mode == "somatic":
            self.mode = mode
        elif mode == "germline":
            self.mode = mode
        else:
            self.mode = False
    def set_bamfile(self, samplename, case="", control=""):
        if case != "":
            self.case = BamFile(case, samplename)
            if not self.case.isexist():
                info(self.case.path + " is not exists! Please check the path!")
        if control!="":
            self.control = BamFile(control, samplename)
            if not self.control.isexist():
                info(self.control.path + " is not exists! Please check the path!")
    def set_config(self, cfg):  #if have two class level, the method can be merge in __init__
        self.cfg = cfg
    def set_seq_type(self, seq_type):
        self.seq_type = seq_type
    def call_variant(self):
        if self.mode == "germline":
            msg = "Running %s in %s mode to %s." % (self.caller_name, self.mode, self.case)
            info (msg)
            out_vcf=""
            if self.caller_name == "tvc":   # For Ion Torrent only
                self.case.index(self.cfg)
                out_vcf = self.case.torrent_caller(self.cfg, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "haplotypecaller":
                self.case.index(self.cfg)
                out_vcf = self.case.haplotype_caller(self.cfg, self.out_dir, seq_type = self.seq_type)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "unifiedgenotyper":
                self.case.index(self.cfg)
                out_vcf = self.case.unifiedgenotyper_caller(self.cfg, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "varscan":
                self.case.index(self.cfg)
                out_vcf = self.case.varscan_caller(self.cfg, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "lofreq":
                self.case.index(self.cfg)
                out_vcf = self.case.lofreq_caller(self.cfg, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "pindel":
                self.case.index(self.cfg)
                out_vcf = self.case.pindel_caller(self.cfg, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if out_vcf:
                msg = "Running %s in %s mode to %s sunccessful!" % (self.caller_name, self.mode, self.case)
                info(msg) 
                return(out_vcf)
            else:
                msg = "Running %s in %s mode to %s fail! Can't generate the %s file." % (
                        self.caller_name, self.mode, self.case, out_vcf)
                return(False)
        elif self.mode == "somatic":
            out_vcf=""
            msg = "Running %s in %s mode to %s and %s." % (self.caller_name, self.mode, self.case, self.control)
            if self.caller_name == "mutect":
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.mutect_caller(self.cfg, self.control, self.out_dir)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "haplotypecaller":
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.haplotype_caller(self.cfg,
                        self.out_dir, self.control, seq_type = self.seq_type)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "unifiedgenotyper":
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.unifiedgenotyper_caller(self.cfg,
                         self.out_dir, self.control.path)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "varscan":
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.varscan_caller(self.cfg,
                         self.out_dir, self.control.path)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "lofreq":
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.lofreq_caller(self.cfg,
                         self.out_dir, self.control.path)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "tvc":   # For Ion Torrent only
                self.case.index(self.cfg)
                self.control.index(self.cfg)
                out_vcf = self.case.torrent_caller(self.cfg, 
                        self.out_dir, self.control.path)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if self.caller_name == "pindel":
                self.case.index(self.samtools)
                out_vcf = self.case.pindel_caller(self.cfg,
                        self.out_dir, self.control.path)
                out_vcf = VcfFile(out_vcf, self.case.samplename)
            if out_vcf:
                msg = "Running %s in %s mode to %s and %s sunccessful!" % (self.caller_name, self.mode, self.case, self.control)
                info(msg)
                return(out_vcf)
            else:
                msg = "Running %s in %s mode to %s fail! Can't generate the %s file." % (
                        self.caller_name, self.mode, self.case, out_vcf)

def variant_caller(options):
    cfg = get_config(options.config)
    caller = cfg["caller"]
    out_dir = options.out_dir
    seq_type = options.seq_type
    samplename = options.samplename
    in_bam = options.in_bam
    if not in_bam:
        return(False)
    caller_list = caller.split(",")
    vcf_list = {}
    for caller in caller_list:
        variantcaller = FundementalCaller(caller)
        variantcaller.set_config(cfg)
        variantcaller.set_seq_type(seq_type)
        variantcaller.set_out_dir(out_dir + "/" + samplename + "/" + caller.capitalize())
        variantcaller.set_bamfile(samplename, in_bam)
        # VcfFile object dict, using all mapped bam file to generate vcf file 
        vcf_list[caller] = variantcaller.call_variant()        
    return(vcf_list)

def variant_caller_somatic(options):
    cfg = get_config(options.config)
    caller = cfg["caller"]
    out_dir = options.out_dir
    seq_type = options.seq_type
    samplename = options.samplename
    casein_bam, controlin_bam = options.in_bam.split(",")
    if not casein_bam or not controlin_bam:
        return(False)
    caller_list = caller.split(",")
    vcf_list = {}
    for caller in caller_list:
        variantcaller = FundementalCaller(caller)
        variantcaller.set_caller_mode("somatic")
        variantcaller.set_config(cfg)
        variantcaller.set_seq_type(seq_type)
        variantcaller.set_out_dir(out_dir + "/" + getid(samplename) + "/" + caller.capitalize())
        variantcaller.set_bamfile(getid(samplename) + "T", casein_bam,controlin_bam)
        # VcfFile object dict, using all mapped bam file to generate vcf file 
        vcf_list[caller] = variantcaller.call_variant()        
    return(vcf_list)


def main():
    options = opt_validate(prepare_optparser())
    if len(options.in_bam.split(",")) ==1:
        variant_caller(options)
    elif len(options.in_bam.split(",")) ==2:
        variant_caller_somatic(options)


if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

