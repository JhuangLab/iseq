import unittest
import os
import tempfile
from iseq.utils import *
from iseq.reffa import *
from iseq.fastq import *
from iseq.csv import *
from iseq.result import *
from iseq.mpileup import *
from iseq.preprocess import *
from iseq.variantcaller import *
from iseq.refinement import *

create_dir("log")
create_dir("restart")
class Tests(unittest.TestCase):
    def setUp(self):
        info('--------- Test setUp Start --------------')
        code_dir = get_root_dir()
        self.code_dir = code_dir
        self.cfg_path = "config.cfg"
        cfg = get_config(self.cfg_path)
        self.cfg = cfg
        self.gatk = cfg["gatk"]
        self.star = cfg["star"]
        self.bwa = cfg["bwa"]
        self.bowtie = cfg["bowtie"] + "-build"
        self.bowtie2 = cfg["bowtie2"] + "-build"
        self.tophat = cfg["tophat"]
        self.reffa = cfg["reffa"]
        self.samtools = cfg["samtools"]
        temp_dir = tempfile.mktemp()
        self.temp_dir = temp_dir
        info("Temp:%s" % temp_dir)
        create_dir(temp_dir)
        info('--------- Test setUp END ----------------')
    def tearDown(self):
        info('--------- Test tearDown Start -----------')
        #rm(self.temp_dir)
        info('--------- Test tearDown END -------------')
    def utils_test(self):
        info('--------- Test utils Module Start -----------')
        temp_dir = self.temp_dir
        #get_config
        self.assertEqual(type(self.cfg), dict)
        self.assertEqual(type(self.star), str)
        self.assertGreater(len(self.star), 3)
        #create_dir & isexist & cp & rm & mv &ln
        self.assertEqual(isexist(temp_dir), True)
        cp(temp_dir, temp_dir + "2")
        self.assertEqual(isexist(temp_dir + "2"), True)
        mv(temp_dir, temp_dir + "3")
        self.assertEqual(isexist(temp_dir), False)
        self.assertEqual(isexist(temp_dir + "3"), True)
        cp(temp_dir + "3", temp_dir)
        ln(temp_dir + "3", temp_dir + "4")
        self.assertEqual(os.path.islink(temp_dir + "4"), True)
        self.assertEqual(rm(temp_dir + "2"), True)
        self.assertEqual(rm(temp_dir + "4"), True)
        self.assertEqual(rm(temp_dir + "3"), True)
        #getid
        self.assertEqual(getid("A04A.fq.gz"), "A04A")
        self.assertEqual(getid("A04A_R1.fq.gz"), "A04A")
        self.assertEqual(getid("A04A_1.fq.gz"), "A04A")
        self.assertEqual(getid("A04A_R2.fq.gz"), "A04A")
        self.assertEqual(getid("A04A_2.fq.gz"), "A04A")
        #get_root_dir
        root_dir = get_root_dir()
        self.assertEqual(type(root_dir), str)

        #gzip and gunzip
        temp_file = temp_dir + "/gunzip.txt"
        fn = open(temp_file, 'a')
        fn.writelines("test")
        fn.close()
        self.assertEqual(isexist(temp_file), True)
        gzip_compress(temp_file, temp_file + ".gz")
        self.assertEqual(isexist(temp_file + ".gz"), True)
        rm(temp_file)
        gzip_uncompress(temp_file + ".gz", temp_file)
        self.assertEqual(isexist(temp_file), True)

        #catmerge
        cp(temp_file, temp_file + "2")
        catmerge(fn_list=[temp_file, temp_file + "2"], out_fn=temp_file + "3" )
        fn = open(temp_file + "3")
        self.assertEqual("".join(fn.readlines()), "testtest")
        fn.close()

        #FundementalFile
        fn_test = FundementalFile(temp_file)
        self.assertEqual(fn_test.basename, os.path.basename(temp_file))
        self.assertEqual(fn_test.runid, "default")
        self.assertEqual(fn_test.dirname, os.path.dirname(temp_file))
        self.assertEqual(fn_test.isexist(), True)
        self.assertEqual(fn_test.cp(temp_file + "3"), True)
        self.assertEqual(isexist(temp_file + "3"), True)
        self.assertEqual(fn_test.mv(temp_file + "2"), True)
        fn_test.path = temp_file + "2"
        self.assertEqual(isexist(temp_file + "2"), True)
        self.assertEqual(isexist(temp_file), False)
        fn_test.ln(temp_file + "4")
        self.assertEqual(isexist(temp_file + "4"), True)
        fn_test.gzip_compress()
        fn_test.path = fn_test.path + ".gz"
        self.assertEqual(isexist(fn_test.path), True)
        rm(temp_file + "2")
        self.assertEqual(fn_test.gzip_uncompress(), True)
        self.assertEqual(isexist(fn_test.path.replace(".gz", "")), True)

        info('--------- Test utils Module END -----------')
    def reffa_test(self):
        info('--------- Test reffa Module Start -----------')
        reffa_path = "example_dat/exampleFASTA.fasta"
        cp(reffa_path, "%s/hg19.fa" % self.temp_dir)
        reffa_path = "%s/hg19.fa" % self.temp_dir
        self.reffa = reffa_path
        pattren = re.compile(".(fa|fasta|FASTA|FA)$")
        replace_str = re.search(pattren, reffa_path).group()
        out_dict = reffa_path.replace(replace_str,".dict")
        reffa_test = ReffaFile(reffa_path, self.cfg)
        reffa_test.generate_dict()
        self.assertEqual(isexist(out_dict), True)
        status = reffa_test.bwa_index()
        self.assertEqual(status, True)
        status = reffa_test.star_index()
        self.assertEqual(status, True)
        status = reffa_test.bowtie_index()
        self.assertEqual(status, True)
        status = reffa_test.bowtie2_index()
        self.assertEqual(status, True)
        status = reffa_test.tmap_index()
        self.assertEqual(status, True)
        info('--------- Test reffa Module END -----------')

    def fastq_test(self):
        info('--------- Test fastq Module Start -----------')
        fastq_1 = "example_dat/reads_1.fq.gz"
        fastq_2 = "example_dat/reads_2.fq.gz"
        cp(fastq_1, "%s/reads_1.fq.gz" % self.temp_dir)
        cp(fastq_2, "%s/reads_2.fq.gz" % self.temp_dir)
        fastq_1 = "%s/reads_1.fq.gz" % self.temp_dir
        fastq_2 = "%s/reads_2.fq.gz" % self.temp_dir
        self.fastq_1 =fastq_1
        self.fastq_2 =fastq_2
        fastq_test_1 = FastqFile(fastq_1, "example", self.cfg)
        fastq_test_2 = FastqFile(fastq_2, "example", self.cfg)
        test_list = ["bwa", "bowtie", "bowtie2", "star", "tophat", "tmap"]
        test_dir = {}
        for i in test_list:
            test_dir_tmp = "%s/%s_test" % (self.temp_dir, i)
            create_dir(test_dir_tmp)
            test_dir.update({i:test_dir_tmp})
        self.cfg["reffa"] = self.reffa
        fastq_test_1.bwa_mapping(test_dir["bwa"], fastq_test_2.path)
        fastq_test_1.bowtie_mapping(test_dir["bowtie"], fastq_test_2.path)
        fastq_test_1.bowtie2_mapping(test_dir["bowtie2"], fastq_test_2.path)
        fastq_test_1.star_mapping(test_dir["star"], fastq_test_2.path)
        fastq_test_1.tophat_mapping(test_dir["tophat"], paired_fastq=fastq_test_2.path)
        fastq_test_1.tmap_mapping(test_dir["tmap"])
        info('--------- Test fastq Module END -----------')

    def bam_test(self):
        info('--------- Test bam Module Start -----------')
        test_bam = "example_dat/exampleBAM.bam"
        cp(test_bam, "%s/exampleBam.bam" % self.temp_dir)
        test_bam = BamFile("%s/exampleBam.bam" % self.temp_dir, "example", config_dict = self.cfg)
        test_bam.index()
        self.cfg["reffa"] = self.reffa
        test_bam.mpileup(test_bam.path + ".mpileup")
        test_bam.contig_reorder(test_bam.path + ".contig_reorder.bam")
        test_bam.add_read_group(test_bam.path + "add_read_group.bam")
        test_bam.mark_duplicates(test_bam.path + "mark.bam", test_bam.path + ".metrics")
        test_bam.realigner_target_creator(test_bam.path + ".intervals")
        test_bam.indel_realigner(test_bam.path + ".intervals", test_bam.path + "inderealgner.bam")
        test_bam.recalibration(test_bam.path + ".grp")
        test_bam.print_reads(test_bam.path + ".grp", test_bam.path + ".final.bam")
        test_grp= "example_dat/exampleBAM.bam.grp"
        test_bam.print_reads(test_grp, test_bam.path + ".final.bam")

        test_bam = "example_dat/exampleBAM.final.bam"
        control_bam = "example_dat/exampleBAM.bam"
        cp(test_bam, "%s/exampleBam.final.bam" % self.temp_dir)
        cp(control_bam, "%s/exampleBam.control.bam" % self.temp_dir)

        test_bam = BamFile("%s/exampleBam.final.bam" % self.temp_dir, "example", config_dict = self.cfg)
        c_bam = BamFile("%s/exampleBam.control.bam" % self.temp_dir, "example", config_dict = self.cfg)
        self.test_bam = test_bam
        self.c_bam = c_bam
        test_bam.index()
        c_bam.index()
        control_bam = c_bam.path
        test_bam.haplotype_caller("%s/test_bam/haplo/" % self.temp_dir)
        test_bam.haplotype_caller("%s/test_bam/haplo/" % self.temp_dir, control_bam)
        test_bam.unifiedgenotyper_caller("%s/test_bam/unified/" % self.temp_dir)
        test_bam.unifiedgenotyper_caller("%s/test_bam/unified/" % self.temp_dir, control_bam)
        test_bam.lofreq_caller("%s/test_bam/lofreq/" % self.temp_dir)
        test_bam.lofreq_caller("%s/test_bam/lofreq/" % self.temp_dir, control_bam)
        test_bam.varscan_caller("%s/test_bam/varscan/" % self.temp_dir)
        test_bam.varscan_caller("%s/test_bam/varscan/" % self.temp_dir, control_bam)
        test_bam.pindel_caller("%s/test_bam/pindel/" % self.temp_dir)
        test_bam.pindel_caller("%s/test_bam/pindel/" % self.temp_dir, control_bam)
        info('--------- Test bam Module END -----------')

    def vcf_test(self):
        test_vcf = "example_dat/exampleBAM.final.vcf"
        control_vcf = "example_dat/exampleBAM.final.control.vcf"
        cp(test_vcf, "%s/exampleBam.final.vcf" % self.temp_dir)
        cp(control_vcf, "%s/exampleBam.final.control.vcf" % self.temp_dir)

        test_vcf = VcfFile("%s/exampleBam.final.vcf" % self.temp_dir, "example")
        control_vcf = VcfFile("%s/exampleBam.final.control.vcf" % self.temp_dir, "example")
        self.test_vcf = test_vcf
        self.control_vcf = control_vcf

        vcf_test_dir = "%s/test_vcf" % self.temp_dir
        create_dir(vcf_test_dir)
        test_vcf.annovar(self.cfg, vcf_test_dir)
        test_vcf.control_filter(self.cfg, control_vcf.path, vcf_test_dir + "/example.filter.vcf")
        test_vcf.control_filter_nosplit(self.cfg, control_vcf.path, vcf_test_dir + "/example.filter.nosplit.vcf")
        test_vcf.fsqd_filter(self.cfg, vcf_test_dir + "example.fsqd.filter.vcf")
        test_vcf.merge(self.cfg, vcf_test_dir + "/example.merge.vcf")
        test_vcf.select_snp(self.cfg, vcf_test_dir + "/example.snp.vcf")
        test_vcf.select_indel(self.cfg, vcf_test_dir + "/example.indel.vcf")
        test_vcf.snpfilter(vcf_test_dir + "/example.snpfilter.vcf")
        test_vcf.unifiedgenotyper_filter(self.cfg, vcf_test_dir + "/example.unifilter.vcf")

    def csv_test(self):
        test_csv = "example_dat/example.avinput.hg19.multianno.csv"
        cp(test_csv, "%s/test.csv" % self.temp_dir)
        test_csv = CsvFile("%s/test.csv" % self.temp_dir, "example")
        test_csv_dir = "%s/test_csv" % self.temp_dir
        create_dir(test_csv_dir)
        test_csv.fmt_annovar(self.cfg, "%s/test_csv_fmt.txt" % test_csv_dir)
        test_csv.get_pos("%s/test_csv.pos" % test_csv_dir)
        test_csv.mpileup(self.cfg, self.test_bam.path, self.test_bam.path + ".mpileup")
        self.mpileup_file = self.test_bam.path + ".mpileup"
        self.result_file = "%s/test_csv_fmt.txt" % test_csv_dir
    
    def mpileup_test(self):
        self.mpileup_file = MpileupFile(self.mpileup_file, "example")
        self.mpileup_file.fmtmpileup(self.mpileup_file.path + ".fmt")
        self.mpileup_file = MpileupFile(self.mpileup_file.path + ".fmt", "example")
        self.snv_frq = self.mpileup_file.path + ".snv.frq"
        self.indel_frq = self.mpileup_file.path + ".indel.frq"
        self.mpileup_file.getSnvFrq(self.snv_frq)
        self.mpileup_file.getIndelFrq(self.indel_frq)
        self.snv_frq_file = MpileupFile(self.snv_frq, "example")
        self.indel_frq_file = MpileupFile(self.indel_frq, "example")

    def result_test(self):
        result_file = ResultFile(self.result_file, "example")
        result_file.fmt_result2final(result_file.path + ".final.txt", self.snv_frq, self.indel_frq)

    def preprocess_test(self):
        cmd = "python %s/iseq/preprocess.py -c %s -s example -1 %s -2 %s -d 1111111111 -o %s" % (self.code_dir, 
                self.cfg_path, self.fastq_1, self.fastq_2, self.temp_dir)
        runcmd(cmd)
        op = iseq_option("example", self.cfg_path, self.fastq_1, self.fastq_2, "1111111", self.test_bam.path, 
                "1111111111", "1", "1", self.temp_dir)
        gp = GenomePreprocessor(op)
        gp.preprocess()
        fq = FastqPreprocessor(op)
        fq.preprocess()
        pre_c = pre_processor(op)
        pre_c.fastq_mapping()
        pre_c.add_read_group()

    def variantcall_test(self):
        op = iseq_option("example", self.cfg_path, self.fastq_1, self.fastq_2, "1111111", self.test_bam.path, 
                "1111111111", "1", "1", self.temp_dir, "dna")
        variant_caller(op)
        op = iseq_option("example", self.cfg_path, self.fastq_1, self.fastq_2, "1111111", 
                self.test_bam.path + "," + self.c_bam.path,
                "1111111111", "1", "1", self.temp_dir, "dna")
        variant_caller_somatic(op)
    
    def refinement_test(self):
        op = iseq_option("example", self.cfg_path, self.fastq_1, self.fastq_2, "1111111", 
                self.test_bam.path,
                "1111111111", "1", "1", self.temp_dir, "dna", self.test_vcf, self.control_vcf, "fsqd_filter,unifiedgenotyper_filter")
        vcf_filter(op)

        op = iseq_option("example", self.cfg_path, self.fastq_1, self.fastq_2, "1111111", 
                self.test_bam.path + "," + self.c_bam.path,
                "1111111111", "1", "1", self.temp_dir, "dna", self.test_vcf, self.control_vcf, "fsqd_filter,unifiedgenotyper_filter.control_filter")
        vcf_filter_somatic(op)


    def test_main(self):
        self.utils_test()
        self.reffa_test()
        self.fastq_test()
        self.bam_test()
        #self.vcf_test()
        #self.csv_test()
        #self.mpileup_test()
        #self.result_test()
        #self.preprocess_test()
        #self.variantcall_test()
        #self.refinement_test()

class iseq_option(object):
    def __init__(self, samplename, config, fastq1, fastq2, bamprocess, 
            inbam, dataprocess, Genomeindex, FastqMapping, outdir, seq_type, case_vcf, control_vcf="", 
            Vcffilter="wespipeline"):
        self.samplename = samplename
        self.config = config
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.bamprocess = bamprocess
        self.in_bam = inbam
        self.case_vcf = case_vcf
        self.control_vcf = control_vcf
        self.dataprocess = dataprocess
        self.Genomeindex = Genomeindex
        self.FastqMapping = FastqMapping
        self.out_dir = outdir
        self.seq_type = seq_type
        self.Vcffilter = Vcffilter
        self.runid = "test"
        self.mode = "fastq2final"

if __name__ == '__main__':
    unittest.main()
