#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@fastq: Fastq File Class
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from bam import *
from sam import *

class FastqFile(FundementalFile):
    """
    Description:
        Class Name: fastq File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
        **_mapping:Use optional bowtie/bowtie2/bwa/tophat-bowtie/tophat-bowtie2/STAR as mapper to conduct mapping step, need point the output file.It has default parmeter,and can be replaced using string.
    """
    def __init__(self, path, samplename, config_dict = "", runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    ################################################# Mapping Step Method ################################################
    def bwa_mapping(self, out_dir, paired_fastq = "", machine = "Hiseq", platform="ILLUMINA", lab="JhuangLab"):
        config_dict = self.config_dict
        mapper = config_dict["bwa"]
        samtools = config_dict["samtools"]
        reffa = config_dict["reffa"]
        thread = config_dict["bwa_thread"]
        SMflag = getid(self.samplename)
        reffa_dir = os.path.dirname(reffa)
        out_bam = out_dir + "/" + self.samplename + ".bam" 
        mapping_runed_bam = BamFile(out_bam ,self.samplename, config_dict) # bwa out_bam is same with mapping_runed_bam
        log = "2> %s/log/%s.Bwa.mapping.log" % (os.getcwd(), self.runid)
        if paired_fastq != "" and isexist(paired_fastq):
            info("Running Bwa mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s mem -t 40 -MR '@RG\\tID:%s-%s\\tPL:%s\\tSM:%s\\tLB:%s' \
                    %s %s %s %s| %s fixmate -O bam - - \
                    | %s sort -@ 30 -T %s -O bam -o %s \
                    " % (mapper, machine, self.samplename, platform, SMflag, \
                          lab, reffa, self.path, paired_fastq, log, samtools, samtools, self.samplename, out_bam)
        else:
            info("Running Bwa mapping step for %s." % (self.path))
            cmd = "%s mem -t 40 -MR '@RG\\tID:%s-%s\\tPL:%s\\tSM:%s\\tLB:%s' \
                    %s %s %s | %s fixmate -O bam - - \
                    | %s sort -@ 30 -T %s -O bam -o %s \
                    " % (mapper, machine, self.samplename, platform, SMflag, \
                          lab, reffa, self.path, log, samtools, samtools, self.samplename, out_bam)
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index()
            return(mapping_runed_bam) 
        else:
            if mapper.lower().find("bwa") >= 0:
                runcmd(cmd)
                savecmd(cmd, self.runid)
                if mapping_runed_bam.isexist():
                    mapping_runed_bam.index()
                    return(mapping_runed_bam) 
                else:
                    return(False)
            else:
                info("Please set correct bwa path!")
    def star_mapping(self, out_dir, paired_fastq = ""):
        config_dict = self.config_dict
        mapper = config_dict["star"]
        reffa = config_dict["reffa"]
        thread = config_dict["star_thread"]
        reffa_dir = os.path.dirname(reffa)
        out_dir = out_dir + "/" 
        out_bam =  out_dir + "Aligned.sortedByCoord.out.bam"
        out_bam = BamFile(out_bam ,self.samplename, config_dict)
        mapping_runed_bam = out_dir + out_bam.samplename + ".bam" 
        mapping_runed_bam = BamFile(mapping_runed_bam ,self.samplename, config_dict) 
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        log = "&> %s/log/%s.Star.mapping.log" % (os.getcwd(), self.runid)
        if paired_fastq != "" and isexist(paired_fastq):
            info("Running STAR mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s --genomeDir %s --readFilesIn %s %s --runThreadN %s --outSAMtype BAM SortedByCoordinate \
                      --outFileNamePrefix %s %s" % (mapper, reffa_dir, self.path, paired_fastq, thread ,out_dir, log) 
        else:
            info("Running STAR mapping step for %s." % (self.path))
            cmd = "%s --genomeDir %s --readFilesIn %s --runThreadN %s --outSAMtype BAM SortedByCoordinate \
                      --outFileNamePrefix %s %s" % (mapper, reffa_dir, self.path, thread ,out_dir, log) 
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index()
            return(mapping_runed_bam) # Exist result bam file before run 
        else:
            if mapper.lower().find("star")>=0:
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct Star path!")
        #Format bam file name to samplename.bam
        if out_bam.isexist():
            if out_bam.mv(mapping_runed_bam.path):
                mapping_runed_bam.index()
                return(mapping_runed_bam) 
            else:
                return(False)
        else:
            return(False)
    def bowtie2_mapping(self, out_dir, paired_fastq = ""):
        config_dict = self.config_dict
        mapper = config_dict["bowtie2"]
        reffa = config_dict["reffa"]
        thread = config_dict["bowtie2_thread"]
        reffa_dir = os.path.dirname(reffa)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        out_sam = out_dir + "/" + self.samplename + ".sam" 
        out_sam = SamFile (out_sam, self.samplename, config_dict)
        out_bam = out_dir + "/" + self.samplename + ".unsorted.bam" # out_bam is formatted file of mapping 
        mapping_runed_bam = out_dir + "/" + self.samplename + ".bam"
        mapping_runed_bam = BamFile(mapping_runed_bam, self.samplename, config_dict)
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        log = "&> %s/log/%s.Bowtie2.mapping.log" % (os.getcwd(), self.runid)
        if paired_fastq != "" and isexist(paired_fastq):
            info("Running Bowtie2 mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s -p %s -x %s -1 %s -2 %s -S %s %s" % (mapper, thread, reffa_prefix, self.path, paired_fastq, out_sam.path, log)
        else:
            info("Running Bowtie2 mapping step for %s." % (self.path))
            cmd = "%s -p %s -x %s -1 %s -S %s %s" % (mapper, thread, reffa_prefix, self.path, out_sam.path, log)
        if not out_sam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("bowtie2"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct bowtie2 path!")
        elif mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index()
            return(mapping_runed_bam) # Exist result bam file before run 
        if out_sam.isexist():
            out_bam = out_sam.convert2bam(out_bam) # run successful the SAM file will be delete
            if out_bam.isexist():
                mapping_runed_bam = out_bam.sort(mapping_runed_bam.path)
                mapping_runed_bam.index()
                out_sam.rm()
                return(mapping_runed_bam) # Not Exist result bam file before run
            else:
                return(False)
        else:
            return(False)
    def bowtie_mapping(self, out_dir, paired_fastq = ""):
        config_dict = self.config_dict
        mapper = config_dict["bowtie"]
        reffa = config_dict["reffa"]
        thread = config_dict["bowtie_thread"]
        reffa_dir = os.path.dirname(reffa)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        out_sam = out_dir + "/" + self.samplename + ".sam" 
        out_sam = SamFile (out_sam, self.samplename, config_dict)
        out_bam = out_dir + "/" + self.samplename + "unsorted.bam" # out_bam is formatted file of mapping 
        mapping_runed_bam = out_dir + "/" + self.samplename + ".bam"
        mapping_runed_bam = BamFile(mapping_runed_bam, self.samplename, config_dict)
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        log = "&> %s/log/%s.Bowtie.mapping.log" % (os.getcwd(), self.runid)
        if paired_fastq != "" and isexist(paired_fastq):
            info("Running Bowtie mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s %s -p %s --chunkmbs 2000  -n 2 -l 28 -e 70 -1 %s -2 %s %s %s" %(mapper, reffa_prefix, thread, self.path, paired_fastq, out_sam.path, log)
        else:
            info("Running Bowtie mapping step for %s." % (self.path))
            cmd = "%s %s -p %s --chunkmbs 2000  -n 2 -l 28 -e 70 -1 %s -2 %s %s %s" %(mapper, reffa_prefix, thread, self.path, out_sam.path, log)
        if not out_sam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("bowtie"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct bowtie path!")
        elif mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index()
            return(mapping_runed_bam) # Exist result bam file before run 
        if out_sam.isexist():
            out_bam = out_sam.convert2bam(out_bam)
            if out_bam.isexist():
                out_bam.sort(mapping_runed_bam.path)
                mapping_runed_bam.index()
                out_sam.rm()
                return(mapping_runed_bam) # Not Exist result bam file before run
            else:
                return(False)
        else:
            return(False)

    def tophat_mapping(self, out_dir, mode = "tophat2", paired_fastq = ""):
        config_dict = self.config_dict
        mapper = config_dict["tophat"]
        reffa = config_dict["reffa"]
        thread = config_dict["tophat_thread"]
        gtf = config_dict["gtf"]
        reffa_dir = os.path.dirname(reffa)
        out_dir = out_dir + "/"
        out_bam =  out_dir + "accepted_hits.bam"
        out_bam = BamFile(out_bam ,self.samplename, config_dict)
        mapping_runed_bam = out_dir + out_bam.samplename + ".bam" # Tophat 1/2 generated out_bam is not same name to mapping_runed_bam file
        mapping_runed_bam = BamFile(mapping_runed_bam, self.samplename, config_dict)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        if mode == "tophat2":
            info("Running tophat2 mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s --num-threads %s --output-dir %s %s %s %s" %(mapper, thread, out_dir, reffa_prefix, self.path, paired_fastq)
            if gtf != "":
                cmd = "%s --num-threads %s -G %s --output-dir %s %s %s %s" %(mapper, thread, gtf, out_dir, reffa_prefix, self.path, paired_fastq)
            log = " &> %s/log/%s.Tophat2.mapping.log" % (os.getcwd(), self.runid)
        else:
            info("Running tophat mapping step for %s and %s." % (self.path, paired_fastq))
            cmd = "%s --num-threads %s --output-dir %s --bowtie1 %s %s %s" %(mapper, thread ,out_dir, reffa_prefix, self.path, paired_fastq)
            if gtf != "":
                cmd = "%s --num-threads %s -G %s --output-dir %s --bowtie1 %s %s %s" %(mapper, thread, gtf, out_dir, reffa_prefix, self.path, paired_fastq)
            log = " &> %s/log/%s.Tophat.mapping.log" % (os.getcwd(), self.runid)
        cmd = cmd + log
        if not out_bam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("tophat"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct tophat path!")
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index() # Exist result bam file before run
            return(mapping_runed_bam)
        if out_bam.isexist():
            if out_bam.mv(out_dir + out_bam.samplename + ".bam"):
                mapping_runed_bam.index() # Not exist result bam file before run
                return(mapping_runed_bam)
            else:
                return(False)
        else:
            return(False)

    def fastqc(self, out_dir):
        config_dict = self.config_dict
        info("Running fastqc for %s." % (self.path))
        fastqc = config_dict["fastqc"]
        path = os.path.expanduser(self.path)
        path = fastq_uncompress(self.path)
        cmd = "%s %s -o %s" % (fastqc, path, out_dir)
        pattren = re.compile(".(fq|fastq)$")
        re_obj = re.search(pattren, path)
        replace_str = re_obj.group()
        filename = os.path.basename(path)
        fastqc_runed_file = "%s/%s" % (out_dir, filename.replace(replace_str, "_fastqc.zip"))
        runcmd(cmd)
        savecmd(cmd, self.runid)
        if isexist(fastqc_runed_file):
            return(True)
        else:
            return(False)
    def __str__(self):
        return(self.path)
def fastq_uncompress(path):
    fn = FundementalFile(path)
    pattren = re.compile(".(GZ|gz)$")
    re_obj = re.search(pattren, path)
    is_gzip = re_obj is not None
    if is_gzip:
        fn.gzip_uncompress()
        replace_str = re_obj.group()
        uncompress_path = path.replace(replace_str, "")
        path = uncompress_path
    return(path)
