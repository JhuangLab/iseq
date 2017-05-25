#! ~/sbin/Python-2.7.7/python
# coding=utf-8
"""Module Description
Copyright (c) 2017 Jianfeng Li <lee_jianfeng@sjtu.edu.cn>
This code is free software; you can redistribute it and/or modify it
under the terms of the MIT License.
@utils Utils module
@status: experimental
@version: 0.0.1
@author: Jianfeng Li
@contact: lee_jianfeng@sjtu.edu.cn
"""
import sys
import string
import re
import math
import os
import time
import subprocess
from optparse import OptionParser
import logging
import time
import copy
import shutil
from pycnf.pycnf import ConfigFile
import gzip
import binascii 
from cStringIO import StringIO
import copy
# filter function
# from filter import filtercreate_dir, indel_poz2gene, snv_poz2gene
# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info
# ------------------------------------
# Config file
# ------------------------------------


def get_config(path="config.cfg"):
    path = os.path.expanduser(path)
    config = ConfigFile(path)
    config.read()
    status = config.merge(config.sections())
    if status:
        return(config.content_merged)
    else:
        return(False)

def create_dir(path, runid = "default"):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        try:
            os.makedirs(path)
            status = True
        except:
            status = False
    else:
        status = True
    if status:
        cmd = "mkdir -p %s" % (path)
        info(cmd)
        savecmd(cmd, runid)
    return(status)

def isexist(fn):
    fn = os.path.expanduser(fn)
    flag = True
    fn = str(fn)
    if not os.path.exists(fn):
        return(False)
    if os.path.getsize(fn) == 0:
        return(False)
    return(True)


def runcmd(cmd, run=True):
    if run:
        info(cmd)
        p = subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        info("Finish runnning command.")
    else:
        info(cmd)


def savecmd(cmd, runid=""):
    if not isexist("restart"):
        os.makedirs("restart")
    fn_name = "restart/%s.sh" % (runid)
    fn = open(fn_name, "a+")
    local_time = time.localtime(time.time())
    format_time = time.strftime('%Y-%m-%d %H:%M',local_time)
    if not format_time in "".join(fn.readlines()):
        cmd = "#" + format_time + "\n" + cmd
    fn.writelines(cmd + "\n")
    fn.close()


def get_root_dir():
    return(os.path.split(os.path.realpath(__file__))[0])


def getid(path):
    path = os.path.expanduser(path)
    path = os.path.basename(path)
    prefix = [".R1.fq.gz", ".R2.fq.gz", ".1_fq.gz", ".2.fq.gz",
              "_R1.fq.gz", "_R2.fq.gz", "_1.fq.gz", "_2.fq.gz"]
    prefix_new = []
    for i in prefix:
        prefix_tmp = i.replace("fq", "fastq")
        prefix_new.append(i)
        prefix_new.append(prefix_tmp)
    for i in [".fq.gz", ".fastq.gz"]:
        prefix_new.append(i)
    sample_id = path
    for j in prefix_new:
        sample_id = sample_id.replace(j, "")
    return(sample_id)

def cp(path, new_path, runid = "default"):
    path = os.path.expanduser(path)
    status = False

    if not isexist(path):
        info(path + " is not exitsts!")
    else:
        is_dir = os.path.isdir(path)
        is_file = os.path.isfile(path)

        if is_file:
            try:
                shutil.copy(path, new_path)
                status = True
            except:
                status = False
        if is_dir:
            try:
                shutil.copytree(path, new_path)
                status = True
            except:
                status = False
        if is_file:
            cmd = "cp %s %s" % (path, new_path)
            info(cmd)
            savecmd(cmd, runid)
        elif is_dir:
            cmd = "cp -r %s %s" % (path, new_path)
            info(cmd)
            savecmd(cmd, runid)
    return(status)

def rm(path, runid = "default"):
    path = os.path.expanduser(path)
    status = False

    if not isexist(path) and ".*" not in path:
        info(path + " is not exitsts!")
    else:
        is_dir = os.path.isdir(path)
        is_file = os.path.isfile(path) or os.path.islink(path)
        if is_file:
            try:
                os.remove(path)
                status = True
            except:
                status = False
        if not is_file and is_dir:
            try:
                shutil.rmtree(path)
                status = True
            except:
                status = False
        if is_file or is_dir:
            cmd = "rm -rf " + path
            info(cmd)
            savecmd(cmd, runid)
    return(status)

def mv(path, new_path, runid = "default"):
    path = os.path.expanduser(path)
    new_path = os.path.expanduser(new_path)
    cmd = "mv %s %s" % (path, new_path)
    if os.path.exists(path):
        os.rename(path, new_path)
        info(cmd)
        savecmd(cmd ,runid) 
        return(True) 
    else:
        info(path + " is not exitsts!")
        return(False)

def ln(path, new_path, runid = "default"):
    path = os.path.expanduser(path)
    new_path = os.path.expanduser(new_path)
    cmd = "ln -s %s %s" % (path, new_path)
    if os.path.exists(path):
        os.symlink(path, new_path)
        info(cmd)
        savecmd(cmd ,runid) 
        return(True) 
    else:
        info(path + " is not exitsts!")
        return(False)

def gzip_compress(raw_data):  
    buf = StringIO()
    f = gzip.GzipFile(mode='wb', fileobj=buf)
    try:
        f.write(raw_data)
    finally:  
        f.close()
    return(buf.getvalue())

def gzip_uncompress(c_data):
    buf = StringIO(c_data)
    f = gzip.GzipFile(mode = 'rb', fileobj = buf)
    try:
        r_data = f.read()
    finally:
        f.close()
    return(r_data)

def __compress_file(fn_in, fn_out):
    f_in = open(fn_in, 'rb')
    f_out = gzip.open(fn_out, 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    if isexist(fn_out):
        return(True)
    else:
        return(False)

def __uncompress_file(fn_in, fn_out):
    f_in = gzip.open(fn_in, 'rb')
    f_out = open(fn_out, 'wb')
    file_content = f_in.read()
    f_out.write(file_content)
    f_out.close()
    f_in.close()
    if isexist(fn_out):
        return(True)
    else:
        return(False)

def gzip_uncompress(path, out_fn = None, runid="default"):
    if out_fn is None:
        out_fn = path.replace(".gz", "")
    cmd = "gunzip -c %s > %s" % (path, out_fn)
    if isexist(out_fn):
        savecmd(cmd, runid) 
        return(True)
    elif isexist(path):
        __uncompress_file(path, out_fn)
        info(cmd)
        savecmd(cmd, runid) 
        return(True)
    else:
        info(path + " is not exitsts!")
        return(False)

def gzip_compress(path, out_fn = None, runid="default"):
    if out_fn is None:
        out_fn = path + ".gz"
    cmd = "gzip -c %s > %s" % (path, out_fn)
    if isexist(out_fn):
        savecmd(cmd, runid) 
        return(True)
    elif isexist(path):
        __compress_file(path, out_fn)
        info(cmd)
        savecmd(cmd, runid) 
        return(True)
    else:
        info(path + " is not exitsts!")

def catmerge(fn_list=[], out_fn = "", runid="default"):
    cmd = "cat"
    out_fn_path = out_fn
    out_fn = open(out_fn, "w")
    dat = ""
    for fn in fn_list:
        fn_tmp = open(fn, "rb")
        dat_tmp = fn_tmp.readlines()
        dat_tmp = "".join(dat_tmp)
        fn_tmp.close()
        dat = dat + dat_tmp 
        cmd = cmd + " " + fn
    out_fn.writelines(dat)
    out_fn.flush()
    if isexist(out_fn_path):
        cmd = cmd + " > " + out_fn_path
        info(cmd)
        savecmd(cmd, runid)
        return(True)
    else:
        return(False)

def set_jdk(config_dict, version = "jdk_17"):
    jdk = config_dict[version]
    os.environ["JAVA_HOME"] = jdk
    os.environ["JRE_HOME"] = jdk + "/jre"
    os.environ["CLASSPATH"] = "%s/lib/dt.jar:%s/lib/tools.jar"
    config_dict["java"] = jdk + "/bin/java"
    return(config_dict)

class FundementalFile(object):
    """
    Description:
        Class Name:All data Class based on this class, need be point the file path to initial.
    Method:
        cp:Copy file or directorie
        rm:Delete the file by cmd "rm path" 
        mv:Rename the file or directorie
        ln:Generate a symlink
        isexist:Check the file wheather exists and file size greater than 0 or not
        gzip_compress:gzip -c path > new.gz 
        gzip_uncompress:gunzip -c path.gz > new
        catmerge:cat file1 file2 ... > newfile
    """
    def __init__(self, path, runid = "default"):
        self.path = str(path)
        self.dirname = os.path.dirname(str(path)) 
        self.basename = os.path.basename (str(path)) 
        self.runid = runid
    def mv(self ,new_path):
        return(mv(self.path, new_path, self.runid))
    def cp(self, new_path):
        return(cp(self.path, new_path, self.runid))
    def rm(self):
        return(rm(self.path, self.runid))
    def ln(self, new_path):
        cmd = "ln -s %s %s" % (self.path, new_path)
        return(ln(self.path, new_path, self.runid))
    def isexist(self):
        return(isexist(self.path))
    def gzip_compress(self, out_fn = None):
        return(gzip_compress(self.path, out_fn, self.runid))
    def gzip_uncompress(self, out_fn = None):
        return(gzip_uncompress(self.path, out_fn, self.runid))
    def catmerge(self, other_fn_list=[], out_fn =""):
        other_fn_list.append(self.path)
        return(catmerge(other_fn_list, out_fn, self.runid))
