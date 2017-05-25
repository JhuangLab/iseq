#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@refinement VCF refinement
@status:  experimental
@version: 0.0.1
@author:  Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from vcf import *
from preprocess import *
from variantcaller import *

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
                         help = "Set the config File.[config.cfg]")
    optparser.add_option("-s", "--samplename", dest = "samplename" ,type = "string",
                         help = "Set the samplename.(Required)")
    optparser.add_option("-f", "--case_vcf", dest = "case_vcf" ,type = "string",
                         help = "You can set this to your case vcf file path.")
    optparser.add_option("-t", "--control_vcf", dest = "control_vcf" ,type = "string",default = "",
                         help = "You can set this to your case vcf file path.")
    optparser.add_option("-i", "--in_bam", dest = "in_bam" ,type = "string",
                         help = "Set your bam file path.(eg. case,control, or case only)")
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
    elif not options.samplename:
        optparser.print_help()
        sys.exit(1)
    elif not options.case_vcf:
        optparser.print_help()
        sys.exit(1)
    return(options)

class FundementalFilter(object):
    def __init__(self, options):
        self.options = options
        self.samplename = options.samplename
        self.case_vcf = VcfFile(options.case_vcf, self.samplename)
        self.out_dir = options.out_dir
        self.in_bam = options.in_bam
        try:
            self.control_vcf = options.control_vcf
        except:
            self.control_vcf = False
        self.cfg = get_config(options.config)
        self.java = self.cfg["java"]
        self.gatk = self.cfg["gatk"]
        self.picard = self.cfg["picard"]
        self.samtools = self.cfg["samtools"]
        self.reffa = self.cfg["reffa"]

class GatkFilter(FundementalFilter):
    """
    Ret:Gatk filter,fsqd_filter is RNAseq specific filter 
    """
    def __init__(self, options):
        FundementalFilter.__init__(self, options)
    def fsqd_filter(self, out_vcf, window = 35, cluster = 3, maxFS = 30.0, minQD = 2.0):  #RNAseq pipeline filter
        out_dir = self.out_dir 
        create_dir(out_dir)
        out_vcf = self.case_vcf.fsqd_filter(self.cfg, out_vcf, window, cluster, maxFS, minQD)
        return(out_vcf)
    def control_filter(self, out_vcf): #wes somatic pipeline unifiedgenotyper_caller drop same of conrol 
        out_dir = self.out_dir
        out_vcf = self.case_vcf.control_filter(self.cfg, self.control_vcf, out_vcf)
        return(out_vcf)
    def unifiedgenotyper_filter(self, out_vcf):
        out_dir = self.out_dir
        out_vcf = self.case_vcf.unifiedgenotyper_filter(self.cfg, out_vcf)
        return(out_vcf)

class Vcffilter(FundementalFilter):
    """
    Ret:Use vcf file filter by gatk, generate annovar file, mpileup file and final standard table
    """
    def __init__(self, options):
        FundementalFilter.__init__(self, options)
        try:
            self.check_filter = options.Vcffilter
        except:
            self.check_filter = "unifiedgenotyper_filter"
        try:
            self.check_annovar = options.vcfannovar
        except:
            self.check_annovar = 1
        try:
            self.check_mpileup = options.mpileup
        except:
            self.check_mpileup = 1
        try:
            self.exononly = options.exononly
        except:
            self.exononly = False
    def gatk_filter(self):   #Rnaseq specific filter in GATK
        pattern = re.compile(".vcf$")
        if "fsqd_filter" in self.check_filter:
            vcffil = GatkFilter(self.options) 
            obj = re.search(pattern, self.case_vcf.path)
            replace_str = obj.group()
            out_vcf =  self.case_vcf.path.replace(replace_str, "_FSQD.vcf")
            self.case_vcf = vcffil.fsqd_filter(out_vcf)
        if "unifiedgenotyper_filter" in self.check_filter: 
            vcffil = GatkFilter(self.options)
            obj = re.search(pattern, self.case_vcf.path)
            replace_str = obj.group()
            out_vcf =  self.case_vcf.path.replace(replace_str, "_UNI.vcf")
            self.case_vcf = vcffil.unifiedgenotyper_filter(out_vcf)
        if "control_filter" in self.check_filter:
            vcffil = GatkFilter(self.options)
            obj = re.search(pattern, self.case_vcf.path)
            replace_str = obj.group()
            out_vcf =  self.case_vcf.path.replace(replace_str, "_CONTR.vcf")
            self.case_vcf = vcffil.control_filter(out_vcf)

    def annovar(self):
        if self.check_annovar:    
            case_vcf = self.case_vcf
            annovar_out_dir = self.case_vcf.dirname + "/annovar/"
            create_dir(annovar_out_dir)
            self.csvfile = case_vcf.annovar(self.cfg, annovar_out_dir)
            out_fmt_annovar_fn = self.csvfile.dirname + "/" + self.samplename + "_result.txt" 
            self.result_fn = self.csvfile.fmt_annovar(self.cfg, out_fmt_annovar_fn)  # Get standard table without mutation frq and depth
    def mpileup(self):
        if self.check_mpileup:
            csvfile = self.csvfile
            out_postion_fn = csvfile.path + ".pos"
            out_mpileup_fn = csvfile.path + ".mpileup"
            out_fmt_mpileup_fn = csvfile.path + ".fmtmpileup"    #Replace blank to unknow
            out_snv_frq_fn = csvfile.path + ".snvfrq"
            out_indel_frq_fn = csvfile.path + ".indelfrq"
            ########### Run Samtools Mpliup #################
            csvfile.get_pos(out_postion_fn, exononly=self.exononly)
            if not isexist(csvfile.pos_file.path):
                info("Annovar or VariantCaller file have 0 item of exon mutation.")
                exit(0)
            mpileup_file = csvfile.mpileup(self.cfg, self.in_bam, out_mpileup_fn)
            mpileup_file = mpileup_file.fmtmpileup(out_fmt_mpileup_fn)
            ########### Get snv and indel alle depth and mutation frquence ############
            self.snv = mpileup_file.get_snv_frq(out_snv_frq_fn)
            self.indel = mpileup_file.get_indel_frq(out_indel_frq_fn)
    def fmt_result2final(self):
        result_fn = self.result_fn
        out_final_fn = self.csvfile.dirname + "/" + self.samplename + ".final.txt"
        self.final_fn = result_fn.fmt_result2final(out_final_fn, self.snv.path, self.indel.path)   # Generate result table
        return(self.final_fn)

class VcffilterSomatic(Vcffilter):
    """
    Ret:Use tumor vcf file filter by gatk, generate annovar file, mpileup file and final standard table
    """
    def mpileup(self):
        casebam = self.in_bam.split(",")[0]
        try:
            controlbam = self.in_bam.split(",")[1]
        except:
            info("please set correct in_bam parameter, eg. --in_bam case.bam,control.bam ")
            sys.exit(1)
        if self.check_mpileup:
            csvfile = self.csvfile
            out_postion_fn = csvfile.path + ".pos"
            out_mpileup_fn = csvfile.path + ".mpileup"
            out_fmt_mpileup_fn = csvfile.path + ".fmtmpileup"    #Replace blank to unknow
            out_snv_frq_fn = csvfile.path + ".snvfrq"
            out_indel_frq_fn = csvfile.path + ".indelfrq"
            ########### Run Samtools Mpliup #################
            csvfile.get_pos(out_postion_fn ,exononly = self.exononly)
            mpileup_file = csvfile.mpileup(self.cfg, casebam, out_mpileup_fn)
            mpileup_file = mpileup_file.fmtmpileup(out_fmt_mpileup_fn)
            ########### Get snv and indel alle depth and mutation frquence ############
            self.casesnv = mpileup_file.get_snv_frq(out_snv_frq_fn)
            self.caseindel = mpileup_file.get_indel_frq(out_indel_frq_fn)
            #Control
            out_mpileup_fn = csvfile.path + ".control.mpileup"
            out_fmt_mpileup_fn = csvfile.path + ".control.fmtmpileup"    #Replace blank to unknow
            out_snv_frq_fn = csvfile.path + ".control.snvfrq"
            out_indel_frq_fn = csvfile.path + ".control.indelfrq"
            mpileup_file = csvfile.mpileup(self.cfg, controlbam, out_mpileup_fn)
            mpileup_file = mpileup_file.fmtmpileup(out_fmt_mpileup_fn)
            self.controlsnv = mpileup_file.get_snv_frq(out_snv_frq_fn)
            self.controlindel = mpileup_file.get_indel_frq(out_indel_frq_fn)
    def fmt_result2final(self):
        result_fn = self.result_fn
        out_final_fn = self.csvfile.dirname + "/" + self.samplename + ".final.txt"
        self.final_fn = result_fn.fmt_result2final(out_final_fn, self.casesnv.path, self.caseindel.path, \
                                                self.controlsnv.path, self.controlindel.path)   # Generate result table
        return(self.final_fn)


#somatic result files filer
def collect_result_file(vcf_out_dir, options, bamfile_list, mode = "somatic", min_case_depth = 10, min_case_frq = 0.08, max_control_frq = 0.05):
    samplename = options.samplename
    root_dir = get_root_dir()
    result_dir = vcf_out_dir+"/result/"
    exon_dir = result_dir + "/exon/"
    config = get_config(options.config)
    caller = config["caller"]
    caller = caller.split(",")
    for m in bamfile_list.keys():
        for c in caller:
            dir_prefix = "%s/%s/%s" % (m.capitalize(), samplename, c.capitalize())
            filename_prefix = "%s.%s.%s" % (samplename, m.capitalize(), c.capitalize())
            old_file_path = "%s/%s/finalResult/%s.txt" % (vcf_out_dir, dir_prefix, 
                    filename_prefix)
            new_file_path = "%s/result/%s" % (vcf_out_dir, filename_prefix)
            create_dir(os.path.dirname(new_file_path))
            cp(old_file_path, new_file_path)
            filename_prefix = "%s.exon.%s.%s" % (samplename, m.capitalize(), c.capitalize())
            old_file_path = "%s/%s/finalResult/%s.txt" % (vcf_out_dir, dir_prefix, 
                filename_prefix)
            new_file_path = "%s/result/exon/%s" %(vcf_out_dir, filename_prefix)
            create_dir(os.path.dirname(new_file_path))
            cp(old_file_path, new_file_path)
    #Filter and collect result file
    exon_input_csv = get_sample_file(exon_dir, samplename)
    cmd = set_filter_cmd(exon_input_csv, exon_dir + samplename + ".exon.collected.txt", mode, 3, 4, min_case_depth, min_case_frq, max_control_frq)
    runcmd(cmd)
    savecmd(cmd, samplename)
    info("The %s all result file can be found in %s!" %(samplename, (vcf_out_dir+"/result")))

def get_sample_file(dirname, samplename):
    files = os.listdir(dirname)
    incsv = ""
    for i in files:
        if (i.split(".")[0]+"x").find(samplename+"x") >= 0 and i.find("collect") < 0:
            if incsv == "":
                incsv = dirname + "/" + i
            else:
                incsv = incsv + "," + dirname + "/" + i
    return(incsv)
def set_filter_cmd(incsv, out_fn , mode = "somatic", mapper_field = 3, caller_field = 4, min_case_depth = 10, min_case_frq = 0.08, max_control_frq = 0.05):
    root_dir = get_root_dir()
    if mode == "somatic":
        cmd = "Rscript %s/Rtools/post_filter.R -i %s --mapperfield %s \
               --callerfield %s --min-casedepth %s --min-casefrq %s \
               --max-controlfrq %s -o %s -m %s" \
                % (root_dir, incsv, mapper_field, caller_field, min_case_depth, \
                    min_case_frq, max_control_frq, out_fn, mode)
    else:
        cmd = "Rscript %s/Rtools/post_filter.R -i %s --mapper_field %s \
               --caller_field %s --min-casedepth %s --min-casefrq %s \
               -o %s -m %s" \
                % (root_dir, incsv, mapper_field, caller_field, min_case_depth, \
                    min_case_frq, out_fn, mode)
    return(cmd)


def vcf_filter(options):
    options.control_vcf = ""
    vcffil = Vcffilter(options)
    vcffil.gatk_filter()
    vcffil.annovar()
    vcffil.mpileup()
    final_fn = vcffil.fmt_result2final()
    return(final_fn)
def vcf_filter_somatic(options):
    vcffil = VcffilterSomatic(options)
    vcffil.gatk_filter()
    vcffil.annovar()
    vcffil.mpileup()
    final_fn = vcffil.fmt_result2final()
    return(final_fn)


def main():
    options = opt_validate(prepare_optparser())
    vcf_filter(options)



if __name__ == "__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0) 

