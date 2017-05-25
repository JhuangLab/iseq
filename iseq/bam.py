#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@bam: Bam File Class
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from reffa import ReffaFile
from vcf import VcfFile
from mpileup import MpileupFile

class BamFile(FundementalFile):
    """
    Description:
        Class Name: BAM File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
        Fundemental Function:
            index:Use samtool index *.bam file and generate *.bam.bai
            mpileup:Use samtool mpileup *.bam file and generate mpileup file
        Preprocesss Function:
            contig_reorder:Use picard to reorder the BAM file according the reference file.
            add_read_group:Use picard to add Read Groups,RGID,RGLB,RGPL,RGPU,RGSM in BAM file header.
            mark_duplicates:Use picard to Mark Duplcates of the BAM file.
            realigner_target_creator:Use GATK to run RealignerTargetCreator
            indel_realigner:Use GATK to run IndelRealigner
            recalibration:Use GATK to run BaseRecalibrator
            print_reads:Use GATK to run PrintReads
            split_ntrim:Use GATK to split_ntrim and conduct ReassignOneMappingQuality.
        Variant Caller:
            haplotype_caller:Use GATK haplotype_caller to conduct Variant Discovery Step.
            unifiedgenotyper_caller:Use GATK unifiedgenotyper_caller to conduct Variant Discovery Step.
            mutect_caller:Use Mutect1 to conduct Variant Discovery Step.
            varscan_caller:Use Varscan2 to conduct Variant Discovery Step.
            torrent_caller:Use TVC to conduct Variant Discovery Step.
            lofreq_caller:Use LoFreq to conduct Variant Discovery Step.
            pindel_caller:Use Pindel to conduct Variant Discovery Step.
    """
    def __init__(self, path, samplename, runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, runid)
        self.samplename = samplename
    ################################################# Fundemental Function ################################################
    def index(self ,config_dict):
        """
        Ret:Use samtool index *.bam file and generate *.bam.bai
        """
        samtools = config_dict["samtools"]
        info("Running the samtools index step for " + self.path)
        cmd = "%s index %s" % (samtools ,self.path)
        if not isexist(self.path + ".bai") and not isexist(self.path[0:-3]+"bai"):
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            savecmd(cmd, self.samplename)
        if isexist(self.path[0:-3]+"bai"):
            info("Run the samtools index step successful!")
            return(self.path[0:-3]+"bai")
        elif isexist(self.path + ".bai"):
            info("Run the samtools index step successful!")
            return(self.path + ".bai")
        else:
            info("Run the samtools index step fail!")
            return(False)

    def mpileup(self, config_dict, out_fn):
        """
        Ret:Use samtool mpileup *.bam file and generate *.bam.mpileup
        """
        samtools = config_dict["samtools"]
        reffa = config_dict["reffa"]
        intervals = config_dict["intervals"]
        info("Running the samtools mpileup step for " + self.path)
        out_fn = MpileupFile(out_fn, self.samplename) 
        if isexist(intervals):
            cmd = "%s mpileup -q 1 -l %s -f %s %s > %s" % (samtools, intervals, reffa, self.path, out_fn.path)
        else:
            cmd = "%s mpileup -q 1 -f %s %s > %s" % (samtools, reffa, self.path, out_fn.path)
        if out_fn.isexist():
            savecmd(cmd, self.samplename)
        elif not self.isexist():
            info("%s BAM file is not exists!" %(self.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd ,self.samplename)
            if not out_fn.isexist():
                info("Samtools mpileup step for %s fail!" %(self.path))
                return(False)
        info("Run the Samtools mpileup step successful!")
        return(out_fn)


    def __str__(self):
        return(self.path)
    ################################################# PreProcesss Function #########################################################
    def contig_reorder(self, config_dict, out_bam):
        """
        Ret:Use picard to reorder the BAM file according the reference file.(BOTH DNA and RNA) 
        """
        reffa = config_dict["reffa"]
        java = config_dict["java"]
        picard = config_dict["picard"]
        info("Now running the contig reorder by picard for %s!" % self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam, self.samplename)
        cmd = "%s -jar %s ReorderSam I=%s O=%s REFERENCE=%s" %(java, picard, in_bam.path, out_bam.path, reffa) 
        if out_bam.isexist():
            savecmd(cmd ,self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd ,self.samplename)
            if not out_bam.isexist():
                info("Contig Reorder step for %s fail!" %(self.path))
                return(False)
        info("Running the contig reorder successful!")
        return(out_bam) # BamFile Class instance 
    def add_read_group(self, config_dict, out_bam, RGID = 1, RGLB = "Jhuanglab", RGPL="ILLUMINA", RGPU = "Hiseq"):
        """
        Ret:Use picard to add Read Groups,RGID,RGLB,RGPL,RGPU,RGSM in BAM file header.(Both DNA and RNA) 
        """
        java = config_dict["java"]
        picard = config_dict["picard"]
        java_max_mem = config_dict["java_max_mem"]
        info("Now running add_read_group step for " + self.path)
        RGSM = self.samplename
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam, self.samplename)
        cmd = "%s -Xmx%s -jar %s AddOrReplaceReadGroups \
                   I=%s \
                   O=%s \
                   RGID=%s \
                   RGLB=%s \
                   RGPL=%s \
                   RGPU=%s \
                   RGSM=%s" %(java, java_max_mem, picard, self.path, out_bam, RGID, RGLB, RGPL, RGPU, RGSM)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("add_read_group step for %s fail!" %(self.path))
                return(False)
        info("Add Read Group run successful!")
        return(out_bam) # BamFile Class instance 
    def mark_duplicates(self, config_dict, out_bam, metrics):
        """
        Ret:Use picard to Mark Duplcates of the BAM file. 
        """
        java = config_dict["java"]
        picard = config_dict["picard"]
        java_max_mem = config_dict["java_max_mem"]
        info("Now running mark_duplicates step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam ,self.samplename) 
        cmd = "%s -Xmx%s -jar %s MarkDuplicates I=%s O=%s  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=%s" %(java, java_max_mem, picard, self.path, out_bam, metrics)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("mark_duplicates step for %s fail!" %(self.path))
                return(False)
        info("mark_duplicates run successful!")
        return(out_bam) # BamFile Class instance 
    def realigner_target_creator(self, config_dict, out_interval):
        """
        Ret:Use GATK realigner_target_creator to the BAM file(DNA Seq use). 
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        thread = config_dict["gatk_preprocess_thread"]
        info("Now running realigner_target_creator step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_interval ,self.samplename)#output is out_interval 
        cmd = "%s -Xmx%s -jar %s  -T RealignerTargetCreator -R %s --num_threads %s \
               -allowPotentiallyMisencodedQuals -I %s \
               -o %s " %(java, java_max_mem, gatk, reffa, thread, self.path, out_interval)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("realigner_target_creator step for %s fail!" %(self.path))
                return(False)
        info("realigner_target_creator run successful!")
        return(out_interval) # BamFile Class instance 
    def indel_realigner(self, config_dict, intervals, out_bam):
        """
        Ret:Use GATK indel_realigner to the BAM file(DNA Seq use). 
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        thread = config_dict["gatk_preprocess_thread"]
        info("Now running indel_realigner step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam ,self.samplename) 
        cmd = "%s -Xmx%s -jar %s  -T IndelRealigner -R %s -targetIntervals %s \
               -allowPotentiallyMisencodedQuals -I %s \
               -o %s " %(java, java_max_mem, gatk, reffa, intervals, self.path, out_bam)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("indel_realigner step for %s fail!" %(self.path))
                return(False)
        info("indel_realigner run successful!")
        return(out_bam) # BamFile Class instance 
    def recalibration(self, config_dict, out_grp):
        """
        Ret:Use GATK recalibration to the BAM file(DNA Seq use). 
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        known_sites_vcf = config_dict["known_sites_vcf"]
        known_sites_vcf = known_sites_vcf.split(":")
        info("Now running recalibration step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_grp ,self.samplename) 
        cmd = "%s -Xmx%s -jar %s  -T BaseRecalibrator -R %s \
               --unsafe -I %s \
               -o %s " %(java, java_max_mem, gatk, reffa, self.path, out_bam)
        for j in known_sites_vcf:
            cmd = cmd + " -knownSites " + j 
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("recalibration step for %s fail!" %(self.path))
                return(False)
        info("recalibration run successful!")
        return(out_grp)
    def print_reads(self, config_dict, grp, out_bam):
        """
        Ret:Use GATK to print DNAseq bam data. 
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        info("Now running print_reads step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam ,self.samplename) 
        cmd = "%s -Xmx%s -jar %s  -T PrintReads -R %s \
               -BQSR %s -I %s \
               -o %s " %(java, java_max_mem, gatk, reffa, grp, self.path, out_bam)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("print_reads step for %s fail!" %(self.path))
                return(False)
        info("print_reads run successful!")
        return(out_bam)

    #split_ntrim: RNA seq bam use
    def split_ntrim(self, config_dict, out_bam):
        """
        Ret:Use GATK to split_ntrim and conduct ReassignOneMappingQuality for RNAseq bam data. 
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        RMQF = config_dict["gatk_RMQF"]
        RMQT = config_dict["gatk_RMQT"]
        java_max_mem = config_dict["java_max_mem"]
        info("Now running splitNtrim step for " + self.path)
        in_bam = BamFile(self.path, self.samplename)
        out_bam = BamFile(out_bam ,self.samplename) 
        cmd = "%s -Xmx%s -jar %s -T SplitNCigarReads \
                  -R %s \
                  -I %s \
                  -o %s \
                  -rf ReassignOneMappingQuality \
                  -RMQF %s \
                  -RMQT %s \
                  -U ALLOW_N_CIGAR_READS" %(java, java_max_mem, gatk, reffa, in_bam, out_bam, RMQF, RMQT)
        if out_bam.isexist():
            savecmd(cmd, self.samplename)
        elif not in_bam.isexist():
            info("%s BAM file is not exists!" %(in_bam.path))
            return(False)
        else:
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_bam.isexist():
                info("split_ntrim step for %s fail!" %(self.path))
                return(False)
        info("split_ntrim run successful!")
        return(out_bam) # BamFile Class instance 

    ######################################################  Variant Caller ###############################################
    def haplotype_caller(self, config_dict, out_dir, control_bam = "", seq_type="dna"):
        """
        Ret:Use GATK HaplotypeCaller to conduct Variant Discovery Step.  
        """
        intervals = config_dict["intervals"]
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        info("Now running haplotype_caller step for " + self.path)
        snp_flag = dbsnp != ""
        intervals_flag = intervals != ""
        create_dir(out_dir)
        out_vcf = out_dir + "/" + self.samplename + ".vcf"
        out_vcf = VcfFile(out_vcf,self.samplename)
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if control_bam != "" and isexist(control_bam):
            if seq_type == "dna":
                cmd = "%s -Xmx%s -Djava.io.tmpdir=%s \
                      -jar %s -R %s \
                      -T HaplotypeCaller \
                      --unsafe\
                      -I %s -I %s -o %s "\
                    % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, control_bam, out_vcf.path)
            else:
                cmd = "%s -Xmx%s -Djava.io.tmpdir=%s \
                      -jar %s -R %s \
                      -T HaplotypeCaller \
                      --unsafe\
                      -dontUseSoftClippedBases \
                      -stand_call_conf 20.0 \
                      -stand_emit_conf 20.0 \
                      -I %s -I %s -o %s "\
                    % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, control_bam, out_vcf.path)

        else:
            if seq_type == "dna":
                cmd = "%s -Xmx%s -Djava.io.tmpdir=%s \
                      -jar %s -R %s \
                      -T HaplotypeCaller \
                      --unsafe\
                      -I %s -o %s"\
                     % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, out_vcf.path)
            else:
                cmd = "%s -Xmx%s -Djava.io.tmpdir=%s \
                      -jar %s -R %s \
                      -T HaplotypeCaller \
                      -dontUseSoftClippedBases \
                      -stand_call_conf 20.0 \
                      -stand_emit_conf 20.0 \
                      --unsafe\
                      -I %s -o %s"\
                     % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, out_vcf.path)

        if self.isexist():
            if not out_vcf.isexist():
                if snp_flag and intervals_flag :
                    cmd = cmd + " --dbsnp %s --intervals %s" %(dbsnp,intervals) 
                elif snp_flag and not intervals_flag:
                    cmd = cmd + " --dbsnp %s" %(dbsnp) 
                elif not snp_flag and intervals_flag:
                    cmd = cmd + " --intervals %s" %(intervals) 
                runcmd(cmd)
                savecmd(cmd, self.samplename)
                if not out_vcf.isexist():
                    return(False)
            else:
                savecmd(cmd , self.samplename)
            info("haplotype_caller caller run successful!")
            return(out_vcf) # VcfFile Class instance
        else:
            info("Bam File not exists, can not conduct haplotype_caller step!")
            return(False)
    def unifiedgenotyper_caller(self, config_dict, out_dir, control_bam = ""):
        """
        Ret:Use GATK UnifiedGenotyper to conduct Variant Discovery Step.  
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        intervals = config_dict["intervals"]
        thread = config_dict["gatk_variantcall_thread"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        create_dir(out_dir)
        def setcmd(bamfile, out_vcf, backrun=False):
            cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -dcov 1000 -nt %s \
                      -T UnifiedGenotyper -glm BOTH \
                      --unsafe \
                      -I %s -o %s "\
                      % (java, java_max_mem, tmp_dir, gatk, reffa, thread, bamfile, out_vcf)
            if snp_flag and intervals_flag :
                cmd = cmd + " --dbsnp %s --intervals %s" %(dbsnp,intervals) 
            elif snp_flag and not intervals_flag:
                cmd = cmd + " --dbsnp %s" %(dbsnp) 
            elif not snp_flag and intervals_flag:
                cmd = cmd + " --intervals %s" %(intervals) 
            if backrun:
                cmd = cmd + " &"
            return(cmd)
        snp_flag = dbsnp != ""
        intervals_flag = intervals != ""
        out_vcf = out_dir + "/" + self.samplename + ".vcf"
        out_vcf = VcfFile(out_vcf, self.samplename)
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if control_bam != "" and isexist(control_bam):
            info("Now running unifiedgenotyper_caller step for " + self.path + " and " + control_bam)
            out_case_vcf = VcfFile(out_vcf.path + ".case", self.samplename)
            out_control_vcf = VcfFile(out_vcf.path + ".control" ,self.samplename)
            case_cmd = setcmd(self.path, out_case_vcf.path)
            control_cmd = setcmd(control_bam, out_control_vcf.path)
            if self.isexist() and isexist(control_bam):
                if not out_vcf.isexist():
                    if not out_case_vcf.isexist():
                        runcmd(case_cmd)
                        savecmd(case_cmd, self.samplename)
                    if not out_control_vcf.isexist():
                        runcmd(control_cmd)
                        savecmd(control_cmd, self.samplename)
                    if not out_case_vcf.isexist() or not out_control_vcf.isexist():
                        info("unifiedgenotyper_caller caller in generate case_vcf or control vcf step fail!")
                        return(False)
                    out_case_vcf.control_filter(config_dict, out_control_vcf.path, out_vcf.path)
                    if not out_vcf.isexist(): 
                        return(False)
                else:
                    savecmd(case_cmd, self.samplename)
                    savecmd(control_cmd, self.samplename)
                    out_case_vcf.control_filter(config_dict, out_control_vcf.path, out_vcf.path)
                    if not out_vcf.isexist(): 
                        return(False)
                info("unifiedgenotyper_caller caller run successful!")
                return(out_vcf) # VcfFile Class instance
            else:
                info("Bam File not exists, can not conduct unifiedgenotyper_caller step!")
                return(False)
        else:
            info("Now running unifiedgenotyper_caller step for " + self.path)
            cmd = setcmd(self.path, out_vcf.path)
            if self.isexist():
                if not out_vcf.isexist():
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                    if not out_vcf.isexist():
                        return(False)
                else:
                    savecmd(cmd, self.samplename)
                info("unifiedgenotyper_caller caller run successful!")
                return(out_vcf) # VcfFile Class instance
            else:
                info("Bam File not exists, can not conduct unifiedgenotyper_caller step!")
                return(False)
    def mutect_caller(self, config_dict, control_bam ,out_dir):
        """
        Ret:Use GATK Mutect to conduct Variant Discovery Step.  
        """
        config_dict = set_jdk(config_dict, "jdk_17")
        java = config_dict["java"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        cosmic = config_dict["cosmic"]
        intervals = config_dict["intervals"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        mutect = config_dict["mutect"]
        create_dir(out_dir)
        info("Now running Mutect step for " + self.path and control_bam)
        snp_flag = dbsnp != ""
        intervals_flag = intervals != ""
        out_vcf = out_dir + "/" + self.samplename + ".vcf"
        tmp = out_vcf + ".tmp"
        out_vcf = VcfFile(out_vcf,self.samplename)
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if isexist(tmp) and not out_vcf.isexist():
            runcmd("grep -v \'REJECT\' %s > %s" % (tmp, out_vcf.path))
        cmd = "%s -jar %s  -T MuTect -R %s -I:tumor %s -I:normal %s \
                --cosmic %s \
                --unsafe \
                -o %s "\
                % (java, mutect, reffa, self.path, control_bam, cosmic, tmp)
        if self.isexist():
            if not out_vcf.isexist() and not isexist(tmp):
                if snp_flag and intervals_flag :
                    cmd = cmd + " --dbsnp %s --intervals %s" %(dbsnp,intervals) 
                elif snp_flag and not intervals_flag:
                    cmd = cmd + " --dbsnp %s" %(dbsnp) 
                elif not snp_flag and intervals_flag:
                    cmd = cmd + " --intervals %s" %(intervals) 
                cmd = cmd + " && grep -v \'REJECT\' %s > %s" % (tmp, out_vcf.path)
                runcmd(cmd)
                savecmd(cmd, self.samplename)
                if not out_vcf.isexist():
                    return(False)
            else:
                savecmd(cmd, self.samplename)
            info("Mutect2 caller run successful!")
            config_dict = set_jdk(config_dict, "jdk_18")
            return(out_vcf) # VcfFile Class instance
        else:
            config_dict = set_jdk(config_dict, "jdk_18")
            info("Bam File not exists, can not conduct mutect_caller step!")
            return(False)
    def varscan_caller(self, config_dict, out_dir="", control_bam = ""):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        varscan = config_dict["varscan"]
        samtools = config_dict["samtools"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        java_max_mem = config_dict["java_max_mem"]
        create_dir(out_dir)
        info("Running varscan_caller step for " + self.path)
        out_snp_vcf = out_dir + "/" + self.samplename + ".snp.vcf"
        out_snp_vcf = VcfFile(out_snp_vcf, self.samplename)
        out_indel_vcf = out_dir + "/" + self.samplename + ".indel.vcf"
        out_indel_vcf = VcfFile(out_indel_vcf, self.samplename)
        out_vcf = out_dir + "/" + self.samplename + ".vcf"
        out_vcf = VcfFile(out_vcf, self.samplename)
        case_bam = BamFile(self.path, self.samplename)
        control_bam = BamFile(control_bam, self.samplename)
        cmd = ""
        if self.isexist():
            if not out_vcf.isexist() and (not out_snp_vcf.isexist() or not out_indel_vcf.isexist()):
                case_mpileup_fn = MpileupFile(out_dir + "/" + self.samplename + ".mpileup.case", self.samplename)
                case_bam.mpileup(config_dict, case_mpileup_fn.path)
                if control_bam.path != "" and control_bam.isexist():
                    control_mpileup_fn = MpileupFile(out_dir + "/" + self.samplename + ".mpileup.control", self.samplename)
                    control_bam.mpileup(config_dict, control_mpileup_fn.path)
                    cmd = "%s -Xmx%s -jar %s somatic %s %s --output-snp %s --output-indel %s --output-vcf"\
                        %(java, java_max_mem, varscan, case_mpileup_fn.path, control_mpileup_fn.path, out_snp_vcf.path, out_indel_vcf.path)
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                else:
                    cmd = "%s -Xmx%s -jar %s mpileup2snp %s --output-vcf 1 --min-var-freq 0.04 > %s"\
                        %(java, java_max_mem, varscan, case_mpileup_fn.path, out_snp_vcf.path)
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                    cmd = "%s -Xmx%s -jar %s mpileup2indel %s --output-vcf 1 --min-var-freq 0.04 > %s"\
                        %(java, java_max_mem, varscan, case_mpileup_fn.path, out_indel_vcf.path)
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                if not out_snp_vcf.isexist() or not out_indel_vcf.isexist():
                    return(False)
                else:
                    out_snp_vcf.varscan2gatkfmt()
                    out_indel_vcf.varscan2gatkfmt()
                    out_snp_vcf.merge(config_dict, out_vcf.path, indel=out_indel_vcf.path)
            else:
                savecmd(cmd, self.samplename)
                out_snp_vcf.varscan2gatkfmt()
                out_indel_vcf.varscan2gatkfmt()
                out_snp_vcf.merge(config_dict, out_vcf.path, indel=out_indel_vcf.path)
            return(out_vcf) # VcfFile Class instance
        else:
            info("Bam File not exists, can not conduct varscan_caller step!")
            return(False)
    def torrent_caller(self, config_dict, out_dir, control_bam=""):
        """
        Ret:Use TVC-5.0.3 to conduct Variant Discovery Step.  
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        tvc = config_dict["tvc"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        intervals = config_dict["intervals"]
        tmp_dir = config_dict["tvc_tmp_dir"]
        thread = config_dict["tvc_thread"]
        json = config_dict["tvc_json"]
        create_dir(out_dir)

        runed_vcf = out_dir + "/" + self.samplename + ".vcf"
        runed_vcf = VcfFile(runed_vcf,self.samplename)
        def setcmd(bamfile, reffa, out_dir, json ="", backrun=False):
            cmd = "%s -i %s -r %s -o %s " \
                      % (tvc, bamfile, reffa, out_dir)
            if json != "":
                cmd = cmd + " -p %s" %(json)
            if backrun:
                cmd = cmd + " &"
            return(cmd)
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if control_bam != "" and isexist(control_bam):
            info("Now running TorrentVariantCaller step for " + self.path + " and " + control_bam)
            out_case_vcf = VcfFile(out_dir + "/case/TSVC_variants.vcf", self.samplename)
            out_control_vcf = VcfFile(out_dir + "/control/TSVC_variants.vcf" ,self.samplename)
            case_cmd = setcmd(self.path, reffa, out_case_vcf.dirname, json)
            control_cmd = setcmd(control_bam, reffa, out_control_vcf.dirname, json)
            if self.isexist() and isexist(control_bam):
                if not runed_vcf.isexist():
                    if not out_case_vcf.isexist():
                        runcmd(case_cmd)
                        savecmd(case_cmd, self.samplename)
                    if not out_control_vcf.isexist():
                        runcmd(control_cmd)
                        savecmd(control_cmd, self.samplename)
                    if not out_case_vcf.isexist() or not out_control_vcf.isexist():
                        info("TorrentVariantCaller caller in generate case_vcf or control vcf step fail!")
                        return(False)
                    out_case_vcf.control_filter(config_dict, out_control_vcf.path, runed_vcf.path)
                    if not runed_vcf.isexist(): 
                        return(False)
                else:
                    savecmd(case_cmd, self.samplename)
                    savecmd(control_cmd, self.samplename)
                    out_case_vcf.control_filter(config_dict, out_control_vcf.path, runed_vcf.path)
                    if not runed_vcf.isexist(): 
                        return(False)
                info("TorrentVariantCaller caller run successful!")
                return(runed_vcf) # VcfFile Class instance
            else:
                info("Bam File not exists, can not conduct TorrentVariantCaller step!")
                return(False)
        else:
            info("Running TorrentVariantCaller step for " + self.path)
            out_vcf= out_dir + "/TSVC_variants.vcf"
            out_vcf = VcfFile(out_vcf, self.samplename) 
            cmd = setcmd(self.path, reffa, out_dir, json)
            if out_vcf.isexist():
                out_vcf.mv(runed_vcf, self.samplename)
            if self.isexist():
                if not runed_vcf.isexist():
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                    if out_vcf.isexist():
                        if not out_vcf.mv(runed_vcf, self.samplename):
                            return(False)
                    else:
                        return(False)
                else:
                    savecmd(cmd, self.samplename)
                info("TorrentVariantCaller run successful!")
                return(runed_vcf) # VcfFile Class instance
            else:
                info("Bam File not exists, can not conduct TorrentVariantCaller step!")
                return(False)
    def lofreq_caller(self, config_dict, out_dir, control_bam = ""):
        """
        Ret:Use lofreq to conduct Variant Discovery Step.  
        """
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        lofreq = config_dict["lofreq"]
        reffa = config_dict["reffa"]
        dbsnp = config_dict["dbsnp"]
        intervals = config_dict["intervals"]
        thread = config_dict["lofreq_thread"]
        create_dir(out_dir)

        info("Running lofreq_caller step for " + self.path)
        out_fn = out_dir + "/" + self.samplename + "_"
        out_snp_vcf = out_dir + "/" + self.samplename + "_somatic_final.snvs.vcf"
        out_indel_vcf = out_dir + "/" + self.samplename + "_somatic_final.indels.vcf"
        runed_vcf = out_dir + "/" + self.samplename + ".vcf"
        runed_vcf = VcfFile(runed_vcf,self.samplename)
        out_snp_vcf = VcfFile(out_snp_vcf, self.samplename) 
        out_indel_vcf = VcfFile(out_indel_vcf, self.samplename) 
        out_snp_vcfgz = FundementalFile(out_snp_vcf.path + ".gz") 
        out_indel_vcfgz = FundementalFile(out_indel_vcf.path + ".gz") 
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if control_bam != "" and isexist(control_bam):
            cmd = "%s somatic -n %s -t %s -f %s -d %s --threads %s --call-indels -o %s" \
                % (lofreq, control_bam, self.path, reffa, dbsnp, thread, out_fn)
            if intervals != "" and isexist(intervals):
                cmd = cmd + " -l %s"%(intervals)
        else:
            cmd = "%s call-parallel --pp-threads %s -f %s --call-indels -o %s " %(lofreq, thread, reffa, runed_vcf)
            if intervals != "" and isexist(intervals):
                cmd = cmd + " -l %s %s"%(intervals, self.path)
            else:
                cmd = cmd + self.path
        if self.isexist():
            if control_bam == "" or (not isexist(control_bam)):
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                if out_snp_vcfgz.isexist() and not out_snp_vcf.isexist():
                    out_snp_vcfgz.gzip_uncompress()
                if out_indel_vcfgz.isexist() and not out_indel_vcf.isexist():
                    out_indel_vcfgz.gzip_uncompress()
                if not runed_vcf.isexist() and out_snp_vcf.isexist() and out_indel_vcf.isexist():
                    out_snp_vcf.merge(config_dict, runed_vcf, indelvcf = out_indel_vcf.path)
                if not runed_vcf.isexist():
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                    out_snp_vcfgz.gzip_uncompress()
                    out_indel_vcfgz.gzip_uncompress()
                    out_snp_vcf.merge(config_dict, runed_vcf, indelvcf = out_indel_vcf.path)
            if runed_vcf.isexist():
                info("Lofreq VariantCaller run successful!")
                return(runed_vcf)
            else:
                return(False)
                info("Lofreq VariantCaller run fail!")
        else:
            info("Bam File not exists, can not conduct lofreq_caller step!")
            return(False)
    def pindel_caller(self, config_dict, out_dir, control_bam=""):
        """
        Ret:Use Pindel to conduct SVs Discovery Step.  
        """
        reffa = config_dict["reffa"]
        pindel_dir = config_dict["pindel_dir"]
        thread = config_dict["pindel_thread"]
        genome_name = config_dict["pindel_genome_name"]
        genome_date = config_dict["pindel_genome_date"]
        insertsize = config_dict["pindel_insertsize"]
        create_dir(out_dir)

        pindel = pindel_dir + "/pindel"
        pindel2vcf4tcga =  pindel_dir + "/pindel2vcf4tcga"
        def __pindelout2vcf(datadir, prefix, out_vcf):
            out_type_list = ["_D","_BP","_SI","_INV","_TD","_LI","_BP"]
            out_fnlist = [ prefix + i for i in out_type_list]
            fn = FundementalFile("/dev/null")
            if not isexist(out_vcf + ".pindelout"):
                fn.catmerge(out_fnlist, out_vcf + ".pindelout")
            cmd = "%s -p %s -r %s -R %s -d %s -v %s -G -so true" \
                    %(pindel2vcf4tcga, out_vcf + ".pindelout", reffa, genome_name, genome_date, out_vcf)
            if not isexist(out_vcf):
                runcmd(cmd)
            savecmd(cmd, self.samplename)
        info("Running Pindel step for " + self.path)
        runed_vcf = VcfFile(out_dir + "/" + self.samplename + ".vcf", self.samplename)
        if isinstance(control_bam, BamFile):
            control_bam = control_bam.path
        if control_bam != "" and isexist(control_bam):
            config_case = out_dir + "/pindel.case.config"
            config_casefn = open(config_case,"w")
            config_casefn.write(self.path + "\t" + insertsize + "\t" + self.samplename + "\n")
            config_casefn.flush()
            config_control = out_dir + "/pindel.control.config"
            config_controlfn = open(config_control,"w")
            config_controlfn.write(control_bam + "\t" + insertsize + "\t" + self.samplename + "\n")
            config_controlfn.flush()
            out_case = out_dir + "/" + self.samplename + ".case"
            out_control = out_dir + "/" + self.samplename + ".control"
            case_cmd = "%s -f %s -i %s -c ALL --number_of_threads %s -o %s" %(pindel, reffa, config_case, thread, out_case) 
            control_cmd = "%s -f %s -i %s -c ALL --number_of_threads %s -o %s" %(pindel, reffa, config_control, thread, out_control) 
        else:
            config_case = out_dir + "/pindel.case.config"
            config_casefn = open(config_case,"w")
            config_casefn.write(self.path + "\t" + insertsize + "\t" + self.samplename + "\n")
            config_casefn.flush()
            out_case = out_dir + "/" + self.samplename + ".case"
            case_cmd = "%s -f %s -i %s -c ALL --number_of_threads %s -o %s" %(pindel, reffa, config_case, thread, out_case) 
        if self.isexist():
            if control_bam != "" and isexist(control_bam):
                if not isexist(out_case + "_D"):
                    runcmd(case_cmd)
                savecmd(case_cmd, self.samplename)
                if not isexist(out_control + "_D"):
                    runcmd(control_cmd)
                savecmd(control_cmd, self.samplename)
                out_case_vcf = VcfFile(out_case + ".vcf", self.samplename)
                out_control_vcf = VcfFile(out_control + ".vcf", self.samplename)
                __pindelout2vcf(out_dir, out_case, out_case_vcf.path)
                __pindelout2vcf(out_dir, out_control, out_control_vcf.path)
                out_case_vcf.control_filter_nosplit(config_dict, out_control_vcf.path, runed_vcf.path)
            else:
                if not isexist(out_case + "_D"):
                    runcmd(case_cmd)
                savecmd(case_cmd, self.samplename)
                out_case_vcf = VcfFile(out_case + ".vcf", self.samplename)
                __pindelout2vcf(out_dir, out_case, out_case_vcf.path)
                out_case_vcf.mv(runed_vcf, self.samplename)  
            if runed_vcf.isexist():
                info("Pindel VariantCaller run successful!")
                return(runed_vcf)
            else:
                return(False)
                info("Pindel VariantCaller run fail!")
        else:
            info("Bam File not exists, can not conduct Pindel step!")
            return(False)
    def __str__(self):
        return(self.path)
