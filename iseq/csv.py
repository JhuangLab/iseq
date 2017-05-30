#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@csv: Csv File Class (ANNOVAR output)
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from result import *
from mpileup import *
class CsvFile(FundementalFile):
    """
    Description:
        Class Name: csv File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
        fmt_annovar: Convert ANNOVAR output file to tab split file, and reduce the colnum
        get_pos: Get the file all positions
        mpileup: Using samtools to mpileup for the all positions
        rbind: Merge CSV file with another CSV file
    """
    def __init__(self, path, samplename, config_dict = "", runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    def fmt_annovar(self, out_fn):
        config_dict = self.config_dict
        root_dir = get_root_dir()
        colnames_need = config_dict["annovar_colnames"]
        out_fn = ResultFile(out_fn, self.samplename, config_dict)
        info("Running format annovar file %s to %s" % (self.path, out_fn.path))
        cmd = "%s %s/Rtools/fmtannovar2result.R -i %s -c %s -o %s -s %s" % (config_dict["Rscript"], root_dir, self.path, colnames_need, out_fn.path, self.samplename)
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_fn.isexist():
                info("Format annovar file %s fail!" % self.path)
                return(False)
        else:
            savecmd(cmd, self.samplename)
        self.fmtfile =  out_fn 
        return(out_fn) # ResultFile Class instance 
    def get_pos(self , out_fn ,sep=",", exononly= False):
        config_dict = self.config_dict
        root_dir = get_root_dir()
        out_fn = MpileupFile(out_fn ,self.samplename, config_dict)
        info("Running get_pos for %s to %s" % (self.path, out_fn.path))
        cmd = "%s %s/Rtools/getPos.R -i %s -o %s --split %s --exononly %s" %(config_dict["Rscript"], root_dir, self.path, out_fn.path, sep, exononly)
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            savecmd(cmd, self.samplename)
        self.pos_file =  out_fn 
        return(out_fn) # MpileupFile Class instance 
    def mpileup(self, in_bam, out_fn):
        config_dict = self.config_dict
        samtools = config_dict["samtools"]
        reffa = config_dict["reffa"]
        root_dir = get_root_dir()
        out_fn = MpileupFile(out_fn ,self.samplename, config_dict)
        info("Running samtools mpileup for %s to %s" % (self.path, out_fn.path))
        cmd = "%s mpileup -q 1 -l %s -f %s %s > %s" % (samtools, self.pos_file ,reffa, in_bam, out_fn.path)
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            savecmd(cmd, self.samplename)
        self.mpileupfile =  out_fn
        return(out_fn) # MpileupFile Class instance 
    def rbind(self, fn, out_fn, paired_header=True, **kwargs):
        fn1 = open(self.path, "r")
        fn2 = open(fn, "r")
        out_fn = open(out_fn,"w")
        for line in fn1:
            out_fn.write(line)
        if paired_header:
            fn2.readline()
        for line in fn2:
            out_fn.write(line)
        for fn in kwargs.values():
            fn2 = open(fn2, "r")
            if paired_header:
                fn2.readline()
            for line in fn2:
                out_fn.write(line)
        if isexist(out_fn):
            return(True)
        else:
            return(False)

    def __str__(self):
        return(self.path)
