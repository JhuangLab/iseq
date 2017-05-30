#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@sam: Sam File Class
@status:  experimental
@version: 0.0.1
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from bam import BamFile

class SamFile(FundementalFile):
    """
    Description:
        Class Name: SAM File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
        conver2bam:Use samtools to convert SAM format file to BAM format file.
    """
    def __init__(self, path, samplename, config_dict = "", runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    def convert2bam(self, out_bam):
        config_dict = self.config_dict
        samtools = config_dict["samtools"]
        info("Running samtools to convert %s SAM File to %s Bam file." % (self.path, out_bam))
        cmd = "%s view -bS %s -o %s" % (samtools, self.path, out_bam)
        out_bam = BamFile(out_bam ,self.samplename, config_dict)
        if  out_bam.isexist():
            savecmd(cmd, self.samplename)
            return(out_bam)
        elif self.isexist():
            if not out_bam.isexist():
                runcmd(cmd)
                savecmd(cmd, self.samplename)
                return(out_bam)
        else:
            info(self.path + " is not exists, cannot conduct the convert2bam step!")
