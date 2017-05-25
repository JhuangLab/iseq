# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@vcf: Vcf File Class
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
from csv import *

class VcfFile(FundementalFile):
    """
    Description:
        Class Name: VCF File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
        fsqd_filter: Use Gatk to filter the VCF file only according the max FS and min QD value
        control_filter_nosplit: Use Gatk SelectVariants to filter the variants existing in the control data (SNP and INDEL no split)
        control_filter: Use Gatk SelectVariants to filter the variants existing in the control data (SNP and INDEL split)
        unifiedgenotyper_filter: A standard filter condition for Unifiedgenotyper output
        snpfilter: Filter the all dbSNP records
        annovar: Use ANNOVAR to annotation the variants
        merge: Use GATK CombineVariants to merge mulitple VCF files
        varscan2gatkfmt: Covert Varscan2 output to GATK format
        select_snp: Select all SNVs mutations
        select_indel: Select all INDELs mutations
    """
    def __init__(self, path, samplename, runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, runid)
        self.samplename = samplename
    def fsqd_filter(self, config_dict, out_vcf, window = 35, cluster = 3, max_FS = 30.0, min_QD = 2.0):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        info("Runing vcf Filter, filter:FS > 30.0, QD < 2.0.")
        out_vcf = VcfFile(out_vcf, self.samplename) 
        cmd = "%s -Xmx%s -jar %s \
                -T VariantFiltration \
                -R %s \
                -V %s \
                -window %s -cluster %s \
                -filterName FS \
                -filter 'FS > %s' \
                -filterName QD \
                -filter 'QD < %s' \
                -o %s " %(java, java_max_mem, gatk, reffa, self.path, window, cluster, max_FS, min_QD, out_vcf.path)
        if self.isexist():
            if not out_vcf.isexist() :
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                savecmd(cmd, self.samplename)
        else:
            info("vcf Filter run fail, " + self.path + " is not exists!")
            return(False)
        info("vcf FS>30 QD<2.0 Filter run successful!")
        return(out_vcf) # VcfFile Class instance
    def control_filter_nosplit(self, config_dict, control_vcf, out_vcf):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        case_vcf = VcfFile(self.path, self.samplename)
        control_vcf = VcfFile(control_vcf, self.samplename)
        out_vcf = VcfFile(out_vcf, self.samplename) 
        def rcmd(case_vcf, control_vcf, out_vcf):
            cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T SelectVariants \
                -V %s --discordance %s -o %s"\
                % (java, java_max_mem, tmp_dir, gatk, reffa, case_vcf , control_vcf, out_vcf)
            if self.isexist():
                if not isexist(out_vcf) :
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                else:
                    savecmd(cmd, self.samplename)
            else:
                info("vcf Filter run fail, " + self.path + " is not exists!")
                return(False)
        if not out_vcf.isexist():
            rcmd(case_vcf.path, control_vcf.path, out_vcf.path)
        info("vcf Control Discordance Filter run successful!")
        return(out_vcf)
    def control_filter(self, config_dict, control_vcf, out_vcf):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        info("Runing vcf Control Filter, select discordance variants compared with control vcf file.")
        case_vcf = VcfFile(self.path, self.samplename)
        control_vcf = VcfFile(control_vcf, self.samplename)
        case_vcf.select_snp(config_dict, self.path + ".snv")
        control_vcf.select_snp(config_dict, control_vcf.path + ".snv")
        case_vcf.select_indel(config_dict, self.path + ".indel")
        control_vcf.select_indel(config_dict, control_vcf.path + ".indel")
        out_vcf = VcfFile(out_vcf, self.samplename) 
        def rcmd(case_vcf_path, control_vcf_path, out_vcf_path):
            cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T SelectVariants \
                -V %s --discordance %s -o %s"\
                % (java, java_max_mem, tmp_dir, gatk, reffa, case_vcf_path , control_vcf_path, out_vcf_path)
            if self.isexist():
                if not isexist(out_vcf_path):
                    runcmd(cmd)
                    savecmd(cmd, self.samplename)
                else:
                    savecmd(cmd, self.samplename)
            else:
                info("vcf Filter run fail, " + self.path + " is not exists!")
                return(False)
        if not out_vcf.isexist():
            if not isexist(case_vcf.path + ".snv.fil.vcf"):
                rcmd(case_vcf.path + ".snv", control_vcf.path + ".snv", case_vcf.path + ".snv.fil.vcf")
            if not isexist(case_vcf.path + ".indel.fil.vcf"):
                rcmd(case_vcf.path + ".indel", control_vcf.path + ".indel", case_vcf.path + ".indel.fil.vcf")
            snvfil = VcfFile(case_vcf.path + ".snv.fil.vcf", self.samplename)
            indelfil = VcfFile(case_vcf.path + ".indel.fil.vcf", self.samplename)
            snvfil.merge(config_dict, out_vcf.path, indelfil = indelfil.path)
            rm(self.path + ".snv*")
            rm(self.path + ".indel*")
            rm(control_vcf.path + ".indel*")
            rm(control_vcf.path + ".snv*")
            rm(out_vcf.path + "*recode.vcf*")
        info("vcf Control Discordance Filter run successful!")
        return(out_vcf) # VcfFile Class instance
    def unifiedgenotyper_filter(self, config_dict, out_vcf):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        info("Runing vcf unifiedgenotyper_filter,add header : 'HARD_TO_VALIDATE,LowCoverage,LowQD,LowQual'.")
        out_vcf = VcfFile(out_vcf, self.samplename) 
        tmp = out_vcf.path + ".tmp"
        cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T VariantFiltration \
              --variant %s \
              --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"  --filterName \"HARD_TO_VALIDATE\" \
              --filterExpression \"DP < 10 \" --filterName \"LowCoverage\" \
              --filterExpression \"QD < 1.5 \" --filterName \"LowQD\" \
              --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \" --filterName \"LowQual\" \
              --clusterWindowSize 10 \
                -o %s \
                && grep '#' %s > %s && grep 'PASS' %s >> %s"\
            % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, tmp, tmp, out_vcf.path, tmp, out_vcf.path)
        if self.isexist():
            if not out_vcf.isexist() :
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                savecmd(cmd, self.samplename)
        else:
            info("vcf Filter run fail, " + self.path + " is not exists!")
            return(False)
        info("vcf unifiedgenotyper_filter run successful!")
        return(out_vcf) # VcfFile Class instance
    def snpfilter(self, out_vcf):
        info("Runing snp filter, it will drop all snp line. ")
        out_vcf = VcfFile(out_vcf, self.samplename) 
        self.cp(out_vcf.path)
        cmd = "sed -i \"/\trs/d\" %s" \
            % (out_vcf.path)
        if self.isexist():
            if not out_vcf.isexist() :
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                savecmd(cmd, self.samplename)
        else:
            info("vcf Filter run fail, " + self.path + " is not exists!")
            return(False)
        info("vcf snp Filter run successful!")
        return(out_vcf) # VcfFile Class instance
    def annovar(self, config_dict, out_dir):
        annovar_dir = config_dict["annovar_dir"]
        buildver = config_dict["annovar_buildver"]
        annovar_flag = config_dict["annovar_flag"]
        info ("Running annovar " + self.path + " VCF file, and output to " + out_dir + ".")
        avinput = out_dir + "/" + self.samplename + ".avinput" 
        out_csv = avinput + "." + buildver + "_multianno.csv"
        out_csv = CsvFile(out_csv, self.samplename)
        create_dir(out_dir)
        cmd = "%s/convert2annovar.pl -withzyg -includeinfo -format vcf4old %s > %s" %(annovar_dir,  self.path, avinput)
        if not isexist(avinput):
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            savecmd(cmd, self.samplename)
        cmd = "%s/table_annovar.pl %s %s/humandb -buildver %s -remove %s -nastring . -csvout" % (annovar_dir, avinput, 
                annovar_dir, buildver, annovar_flag) 
        if not out_csv.isexist() and isexist(avinput):
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        elif out_csv.isexist():
            savecmd(cmd, self.samplename)
            info("Run annovar step for %s successful!" % (self.path))
        else:
            info("Run annovar step for %s fail!" % (self.path))
        return(out_csv) # CsvFile Class instance
    def merge(self, config_dict, out_fn, **kwargs):
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        info ("vcfcombine: Begin to merge vcf outputs for files: %s.", self.path + ", " + ", ".join(kwargs.values()))
        out_fn = VcfFile(out_fn, self.samplename) 
        cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -T CombineVariants \
              -R %s --variant %s -o %s \
              --unsafe --genotypemergeoption UNIQUIFY" \
           % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, out_fn)
        for key in kwargs.keys():
            check = FundementalFile(kwargs[ key ])
            if check.isexist():
                cmd = cmd + " --variant %s " % kwargs[ key ] 
            else:
                info("%s is not exists, it can not be merged.")
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_fn.isexist():
                info("")
                return(False)
        else:
            savecmd(cmd, self.samplename)
        info("Run vcfcombine step for %s successful!" % (self.path + ", " + ", ".join(kwargs.values())))
        return(out_fn) #VcfFile instance
    def varscan2gatkfmt(self):
        root_dir = get_root_dir()
        perl_tools_dir = root_dir + "/Rtools"
        cmd = "perl %s/fmtVarscan2GATK.pl -i %s -o %s" % (perl_tools_dir, self.path, self.path)
        if self.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            info("%s file is not exists!Can not run varscan2gatkfmt step!" % self.path)
            return(False)
        info("Run varscan2gatkfmt step for %s successful!" % (self.path))
    def select_snp(self, config_dict, out_fn):
        vcftools = config_dict["vcftools"]
        runed_vcf = out_fn
        if not isexist (out_fn):
            cmd = "%s --vcf %s --remove-indels --recode --recode-INFO-all --out %s" % (vcftools, self.path, out_fn)
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            runed_vcf = VcfFile(out_fn + ".recode.vcf", self.samplename)
            runed_vcf.mv(out_fn)
        if isexist(out_fn):
            return(True)
        else:
            return(False)
    def select_indel(self, config_dict, out_fn):
        vcftools = config_dict["vcftools"]
        runed_vcf = out_fn
        if not isexist (out_fn):
            cmd = "%s --vcf %s --keep-only-indels --recode --recode-INFO-all --out %s" % (vcftools, self.path, out_fn)
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            runed_vcf = VcfFile(out_fn + ".recode.vcf",self.samplename)
            runed_vcf.cp(out_fn)
        if isexist(out_fn):
            return(True)
        else:
            return(False)
    def __str__(self):
        return(self.path)
