#!/usr/bin/env python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@mpileup: Mpileup File Class
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
from utils import *
class MpileupFile(FundementalFile):
    """
    Description:
        Class Name: mpileup File, need be point the file path and samplename(eg. A101A or A101C and A101T) to initial.
    Method:
    """
    def __init__(self, path, samplename, config_dict = "", runid=None):
        if runid is None:
            runid = samplename
        FundementalFile.__init__(self, path, config_dict, runid)
        self.samplename = samplename
    def fmtmpileup(self, out_fn):
        config_dict = self.config_dict
        rootdir = get_root_dir()
        out_fn = MpileupFile(out_fn, self.samplename, config_dict)
        info("Running fmtmpileup for %s to %s" % (self.path, out_fn.path))
        cmd = "%s %s/Rtools/fmtmpileup.R -i %s -o %s " %(config_dict['Rscript'], rootdir, self.path, out_fn.path)
        if not out_fn.isexist():
            runcmd(cmd)
            savecmd(cmd, self.samplename)
            if not out_fn.isexist():
                info("Formatlize of %s is fail!" % self.path)
                return(False)
        else:
            savecmd(cmd, self.samplename)
        self.fmtmpileupfile =  out_fn
        return(out_fn) # MpileupFile Class instance 
    def get_snv_frq(self, out_fn):
        # Define some dir or file name
        out_fn = MpileupFile(out_fn, self.samplename, self.config_dict) 
        info("Running get_snv_frq for %s to %s" % (self.path, out_fn.path))
        if not out_fn.isexist():
            info ("SNV statistics data will be write in %s" % (out_fn.path))
            fn = open(self.path ,"r")
            fnw = open (out_fn.path,"w")
            snvbase=["A","a","T","t","C","c","G","g","N","n"]
            header="Chr\tPostion\tRef\tmpupDepth\tMpliup\tQuality\tAlleDepth\tCommaNum\tDotNum\tAlt\tAltNum\tFrquency\tAll_Mut:All_MutNum\tAll_MutFrquency\n" 
            fnw.write(header)
            #Calculate counts and max value 
            counts=0
            for line in fn:
                atcgnum = {}
                atcgbothnum = {}
                line=line.strip()
                depth=line.split("\t")[3]
                basecount = self.ComputeBase(line)
                #Get AaTtCcGgNn base count  and both of Aa,Tt,Cc,Gg,Nn from var "basecount"
                count_prefix = "_count"
                for base in snvbase:
                    atcgnum[ base+count_prefix] = basecount[base]
                for base in ["A","T","C","G","N"]:
                    atcgbothnum[base + "_" + base.lower() +count_prefix ] = atcgnum[ base.lower()+ count_prefix ] +atcgnum [ base + count_prefix ] 
                maxitem=""
                for base in atcgbothnum.keys():
                    if maxitem == "":
                        maxitem = base 
                    elif atcgbothnum[ base ] >= atcgbothnum [ maxitem ]:
                        maxitem = base
                altnum = atcgbothnum[maxitem]
                maxitem = maxitem[0]
                #Get other alt and it count
                All_mut=[]
                for base in basecount.keys(): 
                    if base !="comma" or base != "dot":
                        All_mut.append(base) 
                        All_mut.append(str(basecount[base]))
                        All_mut_str=":".join(All_mut)
                #Get snv alle depth
                depth = 0
                for s in All_mut[1::2]:
                    s=str(s)
                    if s.isalnum():
                        depth = int(s) + depth
                All_mut_frq = []
                for num in range(1,len(All_mut),2):
                    if depth == 0:
                        All_mut_frq.append("0")
                    else:
                        All_mut_frq.append(str(float(All_mut[num])/float(depth)))
                    All_mut_frq_str=":".join(All_mut_frq)
                depth = str(depth)
                if depth != "0":
                    row = line + "\t" + str(depth) + "\t" + str(basecount["comma"]) + "\t" + str(basecount["dot"]) +"\t" + maxitem.upper() + "\t" + str(altnum) + "\t"  + str(float(altnum)/float(depth))+ "\t" +All_mut_str + "\t" + All_mut_frq_str + "\n" 
                else:
                    row = line + "\t" + str(depth) + "\t" + str(basecount["comma"]) + "\t" + str(basecount["dot"]) +"\t" + maxitem.upper() + "\t" + str(altnum) +"\t" + "0" + "\t" + All_mut_str + "\t" + All_mut_frq_str + "\n"
                fnw.write(row)
                if(counts%50000==0):
                    info("Total %s row have been processed!" %(str(counts)))
                counts +=1
            fnw.flush()
            info("Total %s row have been processed!" %(str(counts)))
        return(out_fn) # MpileupFile Class instance 

    def get_indel_frq(self,out_fn):
        # Define some dir or file name
        out_fn = MpileupFile(out_fn, self.samplename) 
        info("Running get_indel_frq for %s to %s" % (self.path, out_fn.path))
        if not out_fn.isexist():
            info ("INDEL statistics data will be write in %s" % (out_fn.path))
            fn = open(self.path ,"r")
            fnw = open (out_fn.path,"w")
            snvbase=["A","a","T","t","C","c","G","g","N","n"]
            header="Chr\tPostion\tRef\tmpupDepth\tMpliup\tQuality\tAlleDepth\tCommaNum\tDotNum\tAlt\tAltNum\tFrquency\tAll_Mut:All_MutNum\tAll_MutFrquency\n" 
            fnw.write(header)
            counts=0
            for line in fn:
                atcgnum = {}
                indelnum = {}
                indelbothnum = {}
                line=line.strip()
                basecount = self.ComputeBase(line)
                #Get AaTtCcGgNn base count  and both of Aa,Tt,Cc,Gg,Nn from var "basecount"
                count_prefix = "_count"
                for base in snvbase:
                    atcgnum[ base + count_prefix ] = basecount[ base ]
                #Get INDEL count from var "basecount"
                if len(basecount.keys()) > 12:
                    for base in basecount.keys():
                        if (not base in snvbase) and base != "comma" and base != "dot":
                            indelnum[ base ] = basecount[ base ] 
                    keys = indelnum.keys()
                    for base in keys:
                        try:
                            indelbothnum[ base.upper() ] = indelnum[ base.upper() ] + indelnum[ base.lower() ]
                        except:
                            try:
                                indelbothnum[ base.upper() ] = indelnum[ base.upper() ]
                            except:
                                try:
                                    indelbothnum[ base.upper() ] = indelnum[ base.lower() ]
                                except:
                                    indelbothnum[ base.upper() ] = indelnum[ base]
                #Get Max num in INDEL
                altnum = "0"
                maxitem = "NA"
                if len(basecount.keys()) > 12:
                    for base in indelbothnum.keys():
                        if maxitem == "NA":
                            maxitem = base 
                        elif indelbothnum[ base ] > indelbothnum [ maxitem ]:
                            maxitem = base
                else:
                    altnum = "0" 
                if maxitem == "NA":
                    altnum = "0"
                else:
                    altnum = indelbothnum[maxitem]
                #Get other alt and it count
                All_mut=[]
                for base in basecount.keys(): 
                    if base !="comma" or base != "dot":
                        All_mut.append(base) 
                        All_mut.append(str(basecount[base]))
                        All_mut_str=":".join(All_mut)
                #Get indel alle depth
                depth = 0
                for s in All_mut[1::2]:
                    s=str(s)
                    if s.isalnum():
                        depth = int(s) + depth
                All_mut_frq = []
                for num in range(1,len(All_mut),2):
                    if depth == 0:
                        All_mut_frq.append("0")
                    else:
                        All_mut_frq.append(str(float(All_mut[num])/float(depth)))
                    All_mut_frq_str=":".join(All_mut_frq)
                depth = str(depth)
                if depth != "0":
                    row = line + "\t" + str(depth) + "\t" + str(basecount["comma"]) + "\t" + str(basecount["dot"]) +"\t" + maxitem.upper() + "\t" + str(altnum) + "\t"  + str(float(altnum)/float(depth))+ "\t" +All_mut_str + "\t" + All_mut_frq_str + "\n" 
                else:
                    row = line + "\t" + str(depth) + "\t" + str(basecount["comma"]) + "\t" + str(basecount["dot"]) +"\t" + maxitem.upper() + "\t" + str(altnum) +"\t" + "0" + "\t" + All_mut_str + "\t" + All_mut_frq_str + "\n"
                fnw.write(row)
                if(counts%50000==0):
                    info("Total %s row have been processed in %s!" % (str(counts), out_fn.path))
                counts += 1
            fnw.flush()
            info("Total %s row have been processed!" %(str(counts)))
        return(out_fn) # MpileupFile Class instance 

    #####Process mpileup result file userful function#############
    #Function GetLine: use mpileup result file's one line,return a tuple: (chr,position,reference_alle,mpileup_depth, reads, quality)
    def GetLine(self ,Line="", schr="", sposition=""):
        if isexist(self.path) and len(Line)==0:
            fn=open(self.path,"r")
            for line in fn:
                line=line.strip()
                chrname=line.split("\t")[0]
                position=line.split("\t")[1]
                ref=line.split("\t")[2]
                depth=line.split("\t")[3]
                if depth=="0":
                    reads="-"
                    quality="-"
                else:
                    reads=line.split("\t")[4]
                    quality=line.split("\t")[5]
            if (chrname + position) == (schr + sposition):
                return(chrname,position,ref,depth,reads,quality)
        else:
            chrname=Line.split("\t")[0]
            position=Line.split("\t")[1]
            ref=Line.split("\t")[2]
            depth=Line.split("\t")[3]
            if(depth=="0"):
                reads="-"
                quality="-"
            else:
                reads=Line.split("\t")[4]
                quality=Line.split("\t")[5]
            return(chrname,position,ref,depth,reads,quality)

    # Function ComputeBase: use mpileup result file's one line string, return the all (, . ATCGN INDEL) count.
    def ComputeBase(self, line = "",schr = "",sposition = ""):
        Line = self.GetLine(line ,schr ,sposition)
        allbase=["A","a","T","t","C","c","G","g","N","n"]
        basecount={}
        for i in allbase:
            basecount[str(i)]=0
        commanum = Line[4].count(",")
        dotnum = Line[4].count(".")
        basecount["comma"] = commanum
        basecount["dot"] = dotnum
        INS = re.findall(r"\+[0-9]+[ACGTNacgtn]+", Line[4]) 
        INS_str = " ".join(INS)
        DEL = re.findall(r"-[0-9]+[ACGTNacgtn]+", Line[4]) 
        DEL_str = " ".join(DEL)
        for insbase in list(set(INS)):
            allbase.append(insbase)
        for delbase in list(set(DEL)):
            allbase.append(delbase)
        for base in allbase:
            basecount[base]=Line[4].count(base)
        if(len(allbase)==10):
            return(basecount) 
        else:
            count=0
            for base in allbase:
                if count > 9 :
                   next
                else:
                    basecount[base]=basecount[base] - INS_str.count(base) - DEL_str.count(base)
                    count+=1
        return(basecount) 
    def __str__(self):
        return(self.path)
