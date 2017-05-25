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
    def __init__(self, path, samplename, runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, runid)
        self.samplename = samplename
    ################################################# Mapping Step Method ################################################
    def bwa_mapping(self, config_dict, out_dir, paired_fastq = "", machine = "Hiseq", platform="ILLUMINA", lab="JhuangLab"):
        mapper = config_dict["bwa"]
        samtools = config_dict["samtools"]
        reffa = config_dict["reffa"]
        thread = config_dict["bwa_thread"]
        info("Now running Bwa Mapping step for %s." % (self.samplename))
        SMflag = getid(self.samplename)
        reffa_dir = os.path.dirname(reffa)
        out_bam = out_dir + "/" + self.samplename + ".bam" 
        mapping_runed_bam = BamFile(out_bam ,self.samplename) # bwa out_bam is same with mapping_runed_bam
        if paired_fastq != "" and isexist(paired_fastq):
            cmd = "%s mem -t 40 -MR '@RG\\tID:%s-%s\\tPL:%s\\tSM:%s\\tLB:%s' \
                    %s %s %s | %s fixmate -O bam - - \
                    | %s sort -@ 30 -T %s -O bam -o %s \
                    " % (mapper, machine, self.samplename, platform, SMflag, \
                          lab, reffa, self.path, paired_fastq, samtools, samtools, self.samplename, out_bam)
        else:
            cmd = "%s mem -t 40 -MR '@RG\\tID:%s-%s\\tPL:%s\\tSM:%s\\tLB:%s' \
                    %s %s | %s fixmate -O bam - - \
                    | %s sort -@ 30 -T %s -O bam -o %s \
                    " % (mapper, machine, self.samplename, platform, SMflag, \
                          lab, reffa, self.path, samtools, samtools, self.samplename, out_bam)
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index(config_dict)
            info("Run Bwa Mapping step Successful!")
            return(mapping_runed_bam) 
        else:
            if mapper.lower().find("bwa") >= 0:
                runcmd(cmd)
                savecmd(cmd, self.runid)
                if mapping_runed_bam.isexist():
                    mapping_runed_bam.index(config_dict)
                    info("Run Bwa Mapping step Successful!")
                    return(mapping_runed_bam) 
                else:
                    info("Run Bwa Mapping step fail in generate BAM file!")
                    return(False)
            else:
                info("Please set correct bwa path!")
    def star_mapping(self, config_dict, out_dir, paired_fastq = ""):
        mapper = config_dict["star"]
        reffa = config_dict["reffa"]
        thread = config_dict["star_thread"]
        info("Now running STAR Mapping step for %s." % (self.samplename))
        reffa_dir = os.path.dirname(reffa)
        out_dir = out_dir + "/" 
        out_bam =  out_dir + "Aligned.sortedByCoord.out.bam"
        out_bam = BamFile(out_bam ,self.samplename)
        mapping_runed_bam = out_dir + out_bam.samplename + ".bam" 
        mapping_runed_bam = BamFile(mapping_runed_bam ,self.samplename) 
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        if paired_fastq != "" and isexist(paired_fastq):
            cmd = "%s --genomeDir %s --readFilesIn %s %s --runThreadN %s --outSAMtype BAM SortedByCoordinate \
                      --outFileNamePrefix %s" % (mapper, reffa_dir, self.path, paired_fastq, thread ,out_dir) 
        else:
            cmd = "%s --genomeDir %s --readFilesIn %s --runThreadN %s --outSAMtype BAM SortedByCoordinate \
                      --outFileNamePrefix %s" % (mapper, reffa_dir, self.path, thread ,out_dir) 
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index(config_dict)
            info("Run STAR Mapping step Successful!")
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
                mapping_runed_bam.index(config_dict)
                info("Run STAR Mapping step Successful!")
                return(mapping_runed_bam) 
            else:
                info("Run STAR Mapping step fail in mv BAM file!")
                return(False)
        else:
            info("Run STAR Mapping step fail in generate BAM file!")
            return(False)
    def bowtie2_mapping(self, config_dict, out_dir, paired_fastq = ""):
        mapper = config_dict["bowtie2"]
        reffa = config_dict["reffa"]
        thread = config_dict["bowtie2_thread"]
        info("Now running Bowtie2 Mapping step for %s." % (self.samplename))
        reffa_dir = os.path.dirname(reffa)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        out_sam = out_dir + "/" + self.samplename + ".sam" 
        out_sam = SamFile (out_sam, self.samplename)
        out_bam = out_dir + "/" + self.samplename + ".bam" # out_bam is formatted file of mapping 
        mapping_runed_bam = BamFile(out_bam, self.samplename)
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        if paired_fastq != "" and isexist(paired_fastq):
            cmd = "%s -p %s -x %s -1 %s -2 %s -S %s " % (mapper, thread, reffa_prefix, self.path, paired_fastq, out_sam.path)
        else:
            cmd = "%s -p %s -x %s -1 %s -S %s " % (mapper, thread, reffa_prefix, self.path, out_sam.path)
        if not out_sam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("bowtie2"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct bowtie2 path!")
        elif mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index(config_dict)
            info("Run bowtie_mapping step Successful!")
            return(mapping_runed_bam) # Exist result bam file before run 
        if out_sam.isexist():
            mapping_runed_bam = out_sam.convert2bam(config_dict, out_bam) # run successful the SAM file will be delete
            if mapping_runed_bam.isexist():
                mapping_runed_bam.index(config_dict)
                out_sam.rm()
                info("Run bowtie_mapping step Successful!")
                return(mapping_runed_bam) # Not Exist result bam file before run
            else:
                info("Run bowtie_mapping step fail in generate BAM file!")
                return(False)
        else:
            info("Run bowtie_mapping step fail in generate SAM file!")
            return(False)
    def bowtie_mapping(self, config_dict, out_dir, paired_fastq = ""):
        mapper = config_dict["bowtie"]
        reffa = config_dict["reffa"]
        thread = config_dict["bowtie_thread"]
        info("Now running Bowtie Mapping step for %s." % (self.samplename))
        reffa_dir = os.path.dirname(reffa)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        out_sam = out_dir + "/" + self.samplename + ".sam" 
        out_sam = SamFile (out_sam, self.samplename)
        out_bam = out_dir + "/" + self.samplename + ".bam" # out_bam is formatted file of mapping 
        mapping_runed_bam = BamFile(out_bam, self.samplename)
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        if paired_fastq != "" and isexist(paired_fastq):
            cmd = "%s %s -p %s --chunkmbs 2000  -n 2 -l 28 -e 70 -1 %s -2 %s %s" %(mapper, reffa_prefix, thread, self.path, paired_fastq, out_sam.path)
        else:
            cmd = "%s %s -p %s --chunkmbs 2000  -n 2 -l 28 -e 70 -1 %s -2 %s %s" %(mapper, reffa_prefix, thread, self.path, out_sam.path)
        if not out_sam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("bowtie"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct bowtie path!")
        elif mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index(config_dict)
            info("Run bowtie_mapping step Successful!")
            return(mapping_runed_bam) # Exist result bam file before run 
        if out_sam.isexist():
            mapping_runed_bam = out_sam.convert2bam(config_dict, out_bam)
            if mapping_runed_bam.isexist():
                mapping_runed_bam.index(config_dict)
                out_sam.rm()
                info("Run bowtie_mapping step Successful!")
                return(mapping_runed_bam) # Not Exist result bam file before run
            else:
                info("Run bowtie_mapping step fail in generate BAM file!")
                return(False)
        else:
            info("Run bowtie_mapping step fail in generate SAM file!")
            return(False)

    def tophat_mapping(self, config_dict, out_dir, mode = "tophat2", paired_fastq = ""):
        mapper = config_dict["tophat"]
        reffa = config_dict["reffa"]
        thread = config_dict["tophat_thread"]
        gtf = config_dict["gtf"]
        info("Now running tophat Mapping step for %s." % (self.samplename))
        reffa_dir = os.path.dirname(reffa)
        out_dir = out_dir + "/"
        out_bam =  out_dir + "accepted_hits.bam"
        out_bam = BamFile(out_bam ,self.samplename)
        mapping_runed_bam = out_dir + out_bam.samplename + ".bam" # Tophat 1/2 generated out_bam is not same name to mapping_runed_bam file
        mapping_runed_bam = BamFile(mapping_runed_bam, self.samplename)
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa).group()
        reffa_prefix = reffa.replace(replace_str,"")
        self.path = fastq_uncompress(self.path)
        paired_fastq = fastq_uncompress(paired_fastq)
        if mode == "tophat2":
            cmd = "%s --num-threads %s --output-dir %s %s %s %s" %(mapper, thread, out_dir, reffa_prefix, self.path, paired_fastq)
            if gtf != "":
                cmd = "%s --num-threads %s -G %s --output-dir %s %s %s %s" %(mapper, thread, gtf, out_dir, reffa_prefix, self.path, paired_fastq)
        else:
            cmd = "%s --num-threads %s --output-dir %s --bowtie1 %s %s %s" %(mapper, thread ,out_dir, reffa_prefix, self.path, paired_fastq)
            if gtf != "":
                cmd = "%s --num-threads %s -G %s --output-dir %s --bowtie1 %s %s %s" %(mapper, thread, gtf, out_dir, reffa_prefix, self.path, paired_fastq)
        if not out_bam.isexist() and not mapping_runed_bam.isexist():
            if mapper.lower().find("tophat"):
                runcmd(cmd)
                savecmd(cmd, self.runid)
            else:
                info("Please set correct tophat path!")
        if mapping_runed_bam.isexist():
            savecmd(cmd, self.runid)
            mapping_runed_bam.index(config_dict) # Exist result bam file before run
            info("Run tophat_mapping step Successful!")
            return(mapping_runed_bam)
        if out_bam.isexist():
            if out_bam.mv(out_dir + out_bam.samplename + ".bam"):
                mapping_runed_bam.index(config_dict) # Not exist result bam file before run
                info("Run tophat_mapping step Successful!")
                return(mapping_runed_bam)
            else:
                info("Run tophat_mapping step fail in mv BAM file!")
                return(False)
        else:
            info("Run tophat_mapping step fail in generate BAM file!")
            return(False)

    def fastqc(self, config_dict, out_dir):
        info("Now running fastqc for %s." % (self.samplename))
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
            info("Run fastqc step Successful!")
            return(True)
        else:
            info("Run fastqc step fail!")
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
