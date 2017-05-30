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
        common_filter: A standard filter condition for Unifiedgenotyper output
        snpfilter: Filter the all dbSNP records
        annovar: Use ANNOVAR to annotation the variants
        merge: Use GATK CombineVariants to merge mulitple VCF files
        varscan2gatkfmt: Covert Varscan2 output to GATK format
        select_snp: Select all SNVs mutations
        select_indel: Select all INDELs mutations
    """
    def __init__(self, path, samplename, config_dict = "", runid = None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    def fsqd_filter(self, out_vcf, window = 35, cluster = 3, max_FS = 30.0, min_QD = 2.0):
        config_dict = self.config_dict
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        out_vcf = str(out_vcf)
        info("Runing vcf Filter for %s, filter:FS > 30.0, QD < 2.0." % self.path)
        out_vcf = VcfFile(out_vcf, self.samplename, config_dict) 
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
        log = " &> %s/log/%s.fsqd_filter.log" % (os.getcwd(), self.runid)
        cmd = cmd + log
        if self.isexist():
            if not out_vcf.isexist() :
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                savecmd(cmd, self.samplename)
        else:
            info("vcf Filter run fail, " + self.path + " is not exists!")
            return(False)
        return(out_vcf) # VcfFile Class instance
    def control_filter_nosplit(self, control_vcf, out_vcf):
        config_dict = self.config_dict
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        case_vcf = VcfFile(self.path, self.samplename, config_dict)
        control_vcf = VcfFile(control_vcf, self.samplename, config_dict)
        out_vcf = VcfFile(out_vcf, self.samplename, config_dict) 
        info("Runing vcf Control Filter nosplit for case:%s and control:%s." % self.path, control_vcf)
        def rcmd(case_vcf, control_vcf, out_vcf):
            cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T SelectVariants \
                -V %s --discordance %s -o %s"\
                % (java, java_max_mem, tmp_dir, gatk, reffa, case_vcf , control_vcf, out_vcf)
            log = " &> %s/log/%s.control_filter_nosplit.log" % (os.getcwd(), self.runid)
            cmd = cmd + log
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
        return(out_vcf)
    def control_filter(self, control_vcf, out_vcf):
        config_dict = self.config_dict
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        info("Runing vcf Control Filter split, case:%s and control:%s." % (self.path, control_vcf))
        case_vcf = VcfFile(self.path, self.samplename, config_dict)
        control_vcf = VcfFile(control_vcf, self.samplename, config_dict)
        case_vcf.select_snp(self.path + ".snv")
        control_vcf.select_snp(control_vcf.path + ".snv")
        case_vcf.select_indel(self.path + ".indel")
        control_vcf.select_indel(control_vcf.path + ".indel")
        out_vcf = VcfFile(out_vcf, self.samplename, config_dict) 
        def rcmd(case_vcf_path, control_vcf_path, out_vcf_path, prefix):
            cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T SelectVariants \
                -V %s --discordance %s -o %s"\
                % (java, java_max_mem, tmp_dir, gatk, reffa, case_vcf_path , control_vcf_path, out_vcf_path)
            log = " &> %s/log/%s.%s.control_filter.log" % (os.getcwd(), self.runid, prefix)
            cmd = cmd + log
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
                rcmd(case_vcf.path + ".snv", control_vcf.path + ".snv", case_vcf.path + ".snv.fil.vcf", "snv")
            if not isexist(case_vcf.path + ".indel.fil.vcf"):
                rcmd(case_vcf.path + ".indel", control_vcf.path + ".indel", case_vcf.path + ".indel.fil.vcf", "indel")
            snvfil = VcfFile(case_vcf.path + ".snv.fil.vcf", self.samplename, config_dict)
            indelfil = VcfFile(case_vcf.path + ".indel.fil.vcf", self.samplename, config_dict)
            snvfil.merge(out_vcf.path, indelfil = indelfil.path)
            rm(self.path + ".snv*")
            rm(self.path + ".indel*")
            rm(control_vcf.path + ".indel*")
            rm(control_vcf.path + ".snv*")
            rm(out_vcf.path + "*recode.vcf*")
        return(out_vcf) # VcfFile Class instance
    def common_filter(self, out_vcf):
        config_dict = self.config_dict
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        java_max_mem = config_dict["java_max_mem"]
        info("Running vcf common filter for %s, add header : 'HARD_TO_VALIDATE,LowCoverage,LowQD,LowQual'." % self.path)
        out_vcf = VcfFile(out_vcf, self.samplename, config_dict) 
        tmp = out_vcf.path + ".tmp"
        log = " &> %s/log/%s.Unifiedgenotyper_filter.log" % (os.getcwd(), self.runid)
        cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -R %s -T VariantFiltration \
              --variant %s \
              --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"  --filterName \"HARD_TO_VALIDATE\" \
              --filterExpression \"DP < 10 \" --filterName \"LowCoverage\" \
              --filterExpression \"QD < 1.5 \" --filterName \"LowQD\" \
              --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \" --filterName \"LowQual\" \
              --clusterWindowSize 10 \
                -o %s %s \
                && grep '#' %s > %s && grep 'PASS' %s >> %s"\
            % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, tmp, log, tmp, out_vcf.path, tmp, out_vcf.path)
        if self.isexist():
            if not out_vcf.isexist() :
                runcmd(cmd)
                savecmd(cmd, self.samplename)
            else:
                savecmd(cmd, self.samplename)
        else:
            info("vcf Filter run fail, " + self.path + " is not exists!")
            return(False)
        return(out_vcf) # VcfFile Class instance
    def snpfilter(self, out_vcf):
        config_dict = self.config_dict
        info("Runing snp filter for %s, it will drop all snp line." % self.path)
        out_vcf = VcfFile(out_vcf, self.samplename, config_dict) 
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
        return(out_vcf) # VcfFile Class instance
    def annovar(self, out_dir):
        config_dict = self.config_dict
        annovar_dir = config_dict["annovar_dir"]
        buildver = config_dict["annovar_buildver"]
        annovar_flag = config_dict["annovar_flag"]
        info ("Running annovar " + self.path + " VCF file, and output to " + out_dir + ".")
        avinput = out_dir + "/" + self.samplename + ".avinput" 
        out_csv = avinput + "." + buildver + "_multianno.csv"
        out_csv = CsvFile(out_csv, self.samplename, config_dict)
        create_dir(out_dir)
        cmd = "%s/convert2annovar.pl -withzyg -includeinfo -format vcf4old %s > %s" %(annovar_dir,  self.path, avinput)
        if not isexist(avinput):
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            savecmd(cmd, self.samplename)
        log = " &> %s/log/%s.annovar.log" % (os.getcwd(), self.runid)
        cmd = "%s/table_annovar.pl %s %s/humandb -buildver %s -remove %s -nastring . -csvout" % (annovar_dir, avinput, 
                annovar_dir, buildver, annovar_flag) 
        cmd = cmd + log
        if not out_csv.isexist() and isexist(avinput):
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        elif out_csv.isexist():
            savecmd(cmd, self.samplename)
        else:
            info("Run annovar step for %s fail!" % (self.path))
        return(out_csv) # CsvFile Class instance
    def merge(self, out_fn, **kwargs):
        config_dict = self.config_dict
        java = config_dict["java"]
        gatk = config_dict["gatk"]
        reffa = config_dict["reffa"]
        java_max_mem = config_dict["java_max_mem"]
        tmp_dir = config_dict["gatk_tmp_dir"]
        info ("vcfcombine: Begin to merge vcf outputs for files: %s.", self.path + ", " + ", ".join(kwargs.values()))
        out_fn = VcfFile(out_fn, self.samplename, config_dict) 
        log = " &> %s/log/%s.vcf_merge.log" % (os.getcwd(), self.runid)
        cmd = "%s -Xmx%s -Djava.io.tmpdir=%s -jar %s -T CombineVariants \
              -R %s --variant %s -o %s \
              --unsafe --genotypemergeoption UNIQUIFY" \
           % (java, java_max_mem, tmp_dir, gatk, reffa, self.path, out_fn)
        for key in kwargs.keys():
            check = FundementalFile(kwargs[ key ])
            if check.isexist():
                cmd = cmd + " --variant %s " % kwargs[ key ] 
            else:
                info("%s is not exists, it can not be merged.", check.path)
        cmd = cmd + log
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_fn.isexist():
                info("")
                return(False)
        else:
            savecmd(cmd, self.samplename)
        return(out_fn) #VcfFile instance
    def varscan2gatkfmt(self):
        root_dir = get_root_dir()
        perl_tools_dir = root_dir + "/Rtools"
        info("Running varscan2gatkfmt step for %s." % (self.path))
        cmd = "perl %s/fmtVarscan2GATK.pl -i %s -o %s" % (perl_tools_dir, self.path, self.path)
        if self.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
        else:
            info("%s file is not exists!Can not run varscan2gatkfmt step!" % self.path)
            return(False)
    def select_snp(self, out_fn):
        config_dict = self.config_dict
        vcftools = config_dict["vcftools"]
        runed_vcf = out_fn
        info("Runing select snvs step for %s." % (self.path))
        if not isexist (out_fn):
            cmd = "%s --vcf %s --remove-indels --recode --recode-INFO-all --out %s" % (vcftools, self.path, out_fn)
            log = " &> %s/log/%s.select_snp.log" % (os.getcwd(), self.runid)
            cmd = cmd + log
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            runed_vcf = VcfFile(out_fn + ".recode.vcf", self.samplename, config_dict)
            runed_vcf.mv(out_fn)
        if isexist(out_fn):
            return(True)
        else:
            return(False)
    def select_indel(self, out_fn):
        config_dict = self.config_dict
        vcftools = config_dict["vcftools"]
        runed_vcf = out_fn
        info("Runing select indels step for %s." % (self.path))
        if not isexist (out_fn):
            cmd = "%s --vcf %s --keep-only-indels --recode --recode-INFO-all --out %s" % (vcftools, self.path, out_fn)
            log = " &> %s/log/%s.select_indel.log" % (os.getcwd(), self.runid)
            cmd = cmd + log
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            runed_vcf = VcfFile(out_fn + ".recode.vcf",self.samplename, config_dict)
            runed_vcf.cp(out_fn)
        if isexist(out_fn):
            return(True)
        else:
            return(False)
    def __str__(self):
        return(self.path)
