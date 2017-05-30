#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@result: Result File Class
@status:  experimental
@version: 0.0.1
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *

class ResultFile(FundementalFile):
    """
    Description:
        Class Name: Result txt File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
    """
    def __init__(self, path, samplename, config_dict = "", runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    def fmt_result2final(self, out_fn, case_snp_frq, case_indel_frq, control_snp_frq="", control_indel_frq=""):
        rootdir = get_root_dir()
        out_fn = ResultFile(out_fn, self.samplename, self.config_dict)
        info("Running fmt_result2final for %s to %s." % (self.path, out_fn.path))
        if control_snp_frq == "":
            cmd = "%s %s/Rtools/fmtresult2final.R -i %s --casesnv %s --caseindel %s --output %s" %(self.config_dict["Rscript"], rootdir, self.path, case_snp_frq, case_indel_frq, out_fn.path) 
        else:
            cmd = "%s %s/Rtools/fmtresult2final.R -i %s --casesnv %s --caseindel %s \
                        --controlsnv %s --controlindel %s \
                        --output %s" %(self.config_dict["Rscript"], rootdir, self.path, case_snp_frq, case_indel_frq, control_snp_frq, control_indel_frq, out_fn.path) 

        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_fn.isexist(): 
                info("Run fmt_result2final for %s fail!" % self.path)
                return(False)
        else:
            savecmd(cmd, self.samplename)
        return(out_fn) # ResultFile Class instance
    def __str__(self):
        return(self.path)
