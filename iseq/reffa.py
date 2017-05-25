#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@reffa: Genome Reference File Class
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *

class ReffaFile(FundementalFile):
    """
    Description:
        Class Name: Refference sequence File, need be point the file path to initial.
    Method:
        generate_dict:Use picard to generate genome dict File.
        **_index:Use bwa/bowtie/bowtie2/STAR to generate index of genome File.
    """
    def __init__(self, path, runid = "default"):
        FundementalFile.__init__(self, path, runid)
    def generate_dict(self, config_dict, extra_cmd=""):
        java = config_dict["java"]
        picard = config_dict["picard"]
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        out_dict = self.path.replace(replace_str,".dict")
        cmd = "%s -jar %s CreateSequenceDictionary R=%s O=%s %s" %(java, picard, self.path, out_dict, extra_cmd)
        if not isexist(os.path.expanduser(out_dict)):
            runcmd(cmd)
        savecmd(cmd, self.runid)
        info("Generate dictionary for %s successful!" % self.path)
        if isexist(out_dict):
            return(True)
        else:
            return(False)
    def bwa_index(self, config_dict):
        indexer = config_dict["bwa"]
        thread = config_dict["bwa_thread"]
        info("Now running bwa_index step for " + self.path)
        index_runed_file = self.path + ".bwt" # Set Bwa runned outfile
        cmd = "%s index -a bwtsw %s" %(indexer, self.path)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bwa") >= 0:
                runcmd(cmd)
            else:
                info("Please set correct Bwa path!")
        savecmd(cmd, self.runid)
        info("bwa_index step run successful!")
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def star_index(self, config_dict, mode="1pass", out_dir=None, samplename="", sjtab=""):
        indexer = config_dict["star"]
        thread = config_dict["star_thread"]
        info("Running star_index " + mode + " step for " + self.path)
        if mode == "1pass":
            if out_dir is None:
                out_dir = self.dirname + "/"
            index_runed_file = out_dir + "/SA"
        else:
            if out_dir is None:
                out_dir = self.dirname + "/2pass/" + samplename + "/"
            create_dir(out_dir)
            index_runed_file = out_dir + "/SA" #Set STAR Index Runed File for 1pass or 2pass 
        if mode == "1pass":
            cmd = "%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s  --runThreadN %s" %(indexer, out_dir, self.path, thread)
        else:
            cmd = "%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s  --sjdbFileChrStartEnd %s  --sjdbOverhang 75 --runThreadN %s" %(indexer, out_dir, self.path, sjtab, thread)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("star"):
                runcmd(cmd)
            else:
                info("Please set correct STAR path!")
        savecmd(cmd, self.runid)
        info("STAR index " + mode + " step run successful!")
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def bowtie_index(self, config_dict):
        indexer = config_dict["bowtie"]
        info("Now running bowtie_index step for " + self.path)
        out_dir = self.dirname
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        genome_prefix = self.path.replace(replace_str,"")
        index_runed_file = genome_prefix + ".1.ebwt"
        cmd = "%s %s %s" %(indexer, self.path, genome_prefix)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bowtie"):
                indexer = indexer + "-build"
                cmd = "%s %s %s" %(indexer, self.path, genome_prefix)
                runcmd(cmd)
            else:
                info("Please set correct bowtie path!")
        savecmd(cmd, self.runid)
        info("bowtie_index step run successful!")
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def bowtie2_index(self, config_dict):
        indexer = config_dict["bowtie2"]
        info("Now running bowtie2_index step for " + self.path)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        genome_prefix = self.path.replace(replace_str,"")
        index_runed_file = genome_prefix + ".1.bt2"
        cmd = "%s %s %s" %(indexer, self.path, genome_prefix)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bowtie2"):
                indexer = indexer + "-build"
                cmd = "%s %s %s" %(indexer, self.path, genome_prefix)
                runcmd(cmd)
            else:
                info("Please set correct bowtie2 path!")
        savecmd(cmd, self.runid)
        info("bowtie2_index step run successful!")
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def __str__(self):
        return(self.path)
