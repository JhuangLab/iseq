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
    def __init__(self, path, config_dict = "", runid = "default"):
        FundementalFile.__init__(self, path, config_dict, runid)
    def generate_dict(self, extra_option=None):
        info("Running CreateSequenceDictionary step for " + self.path)
        config_dict = self.config_dict
        java = config_dict["java"]
        picard = config_dict["picard"]
        if extra_option is None:
            extra_option = config_dict["reffafile_generate_dict_extra"]
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        out_dict = self.path.replace(replace_str,".dict")
        log = " &> %s/log/%s.generate_dict.log" % (os.getcwd(), self.runid)
        cmd = "%s -jar %s CreateSequenceDictionary R=%s O=%s %s %s" %(java, picard, self.path, out_dict, extra_option, log)
        if not isexist(os.path.expanduser(out_dict)):
            runcmd(cmd)
        savecmd(cmd, self.runid)
        if isexist(out_dict):
            return(True)
        else:
            return(False)
    def bwa_index(self):
        config_dict = self.config_dict
        indexer = config_dict["bwa"]
        extra_option = config_dict["reffafile_bwa_index_extra"]
        thread = config_dict["reffafile_bwa_index_thread"]
        info("Running Bwa index step for " + self.path)
        index_runed_file = self.path + ".bwt" # Set Bwa runned outfile
        log = "&> %s/log/%s.bwa_index.log" % (os.getcwd(), self.runid)
        cmd = "%s index %s -a bwtsw %s %s" %(indexer, extra_option, self.path, log)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bwa") >= 0:
                runcmd(cmd)
            else:
                info("Please set correct Bwa path!")
        savecmd(cmd, self.runid)
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def star_index(self, mode="1pass", out_dir=None, samplename="", sjtab=""):
        config_dict = self.config_dict
        indexer = config_dict["star"]
        thread = config_dict["reffafile_star_index_thread"]
        pass1_extra_option = config_dict["reffafile_star_index_pass1_extra"]
        pass2_extra_option = config_dict["reffafile_star_index_pass2_extra"]
        info("Running STAR index " + mode + " step for " + self.path)
        if mode == "1pass":
            if out_dir is None:
                out_dir = self.dirname + "/"
            index_runed_file = out_dir + "/SA"
        else:
            if out_dir is None:
                out_dir = self.dirname + "/2pass/" + samplename + "/"
            create_dir(out_dir)
            index_runed_file = out_dir + "/SA" #Set STAR Index Runed File for 1pass or 2pass 
        log = "&> %s/log/%s.Star.index.%s.log" % (os.getcwd(), self.runid, mode)
        if mode == "1pass":
            cmd = "%s %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s  --runThreadN %s %s" %(indexer, pass1_extra_option, out_dir, self.path, thread, log)
        else:
            cmd = "%s %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s  --sjdbFileChrStartEnd %s  --sjdbOverhang 75 --runThreadN %s %s" %(indexer, pass2_extra_option, out_dir, self.path, sjtab, thread, log)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("star"):
                runcmd(cmd)
            else:
                info("Please set correct STAR path!")
        savecmd(cmd, self.runid)
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def bowtie_index(self):
        config_dict = self.config_dict
        indexer = config_dict["bowtie"]
        extra_option = config_dict["reffafile_bowtie_index_extra"]
        info("Running Bowtie index step for " + self.path)
        out_dir = self.dirname
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        genome_prefix = self.path.replace(replace_str,"")
        index_runed_file = genome_prefix + ".1.ebwt"
        log = "&> %s/log/%s.bowtie_index.log" % (os.getcwd(), self.runid)
        cmd = "%s %s %s %s %s" %(indexer, extra_option, self.path, genome_prefix, log)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bowtie"):
                indexer = indexer + "-build"
                cmd = "%s %s %s %s" %(indexer, self.path, genome_prefix, log)
                runcmd(cmd)
            else:
                info("Please set correct bowtie path!")
        savecmd(cmd, self.runid)
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def bowtie2_index(self):
        config_dict = self.config_dict
        indexer = config_dict["bowtie2"]
        extra_option = config_dict["reffafile_bowtie2_index_extra"]
        info("Running Bowtie2 index step for " + self.path)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, self.path).group()
        genome_prefix = self.path.replace(replace_str,"")
        index_runed_file = genome_prefix + ".1.bt2"
        log = "&> %s/log/%s.bowtie2_index.log" % (os.getcwd(), self.runid)
        cmd = "%s %s %s %s %s" %(indexer, extra_option, self.path, genome_prefix, log)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("bowtie2"):
                indexer = indexer + "-build"
                cmd = "%s %s %s %s" %(indexer, self.path, genome_prefix, log)
                runcmd(cmd)
            else:
                info("Please set correct bowtie2 path!")
        savecmd(cmd, self.runid)
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def tmap_index(self):
        config_dict = self.config_dict
        indexer = config_dict["tmap"]
        extra_option = config_dict["reffafile_tmap_index_extra"]
        info("Running TMAP index step for " + self.path)
        index_runed_file = self.path + ".tmap.anno" # Set tmap runned outfile
        log = "&> %s/log/%s.tmap_index.log" % (os.getcwd(), self.runid)
        cmd = "%s index -f %s %s %s" %(indexer, self.path, extra_option, log)
        print(cmd)
        if not isexist(os.path.expanduser(index_runed_file)):
            if indexer.lower().find("tmap") >= 0:
                runcmd(cmd)
            else:
                info("Please set correct tmap path!")
        savecmd(cmd, self.runid)
        if isexist(index_runed_file):
            return(True)
        else:
            return(False)
    def __str__(self):
        return(self.path)
