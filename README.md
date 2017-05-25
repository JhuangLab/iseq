#iseq:A convenient integreted tool of NGS Variant Discovery(RNAseq & DNAseq)
This is the official development repository for iseq.<br>
##Requirements
*JAVA
*Git 2.5 or greater
*GATK3.5
*TVC5.0.2
*Bwa
*STAR_2.4.2a
*Bowtie1.1.2
*Bowtie2-2.2.6
*Tophat2.1.0
*samtools1.3
*Picard1.138
*Annovar
##Install iseq
```Bash
git clone git@bioinfo.rjh.com.cn:lab805/iseq.git #bash
```
##Usage
###RNAseq pipeline
####Fastq files to final result
```Bash
python rnaseq.py -c config.cfg -1 sample_R1.fq -2 sample_R2.fq -m Fastq2vcf -o outdir 
```
####Fastq files to final result in somatic.
```Bash
python rnaseq_somatic.py -c config.cfg -1 case_R1.fq -2 case_R2.fq -3 control_R1.fq -4 control_R2.fq -m Fastq2vcf -o outdir 
```
####Bam files to final result.
```Bash
python rnaseq.py -c config.cfg --inbam case.bam -m Bam2vcf -o outdir 
```
####Bam files to final result in somatic.
```Bash
python rnaseq_somatic.py -c config.cfg --inbam case.bam,control.bam -m Bam2vcf -o outdir 
```
###DNAseq pipeline
####Fastq files to final result
```Bash
python weseq.py -c config.cfg -1 sample_R1.fq -2 sample_R2.fq -m Fastq2vcf -o outdir 
```
####Fastq files to final result in somatic.
```Bash
python weseq_somatic.py -c config.cfg -1 case_R1.fq -2 case_R2.fq -3 control_R1.fq -4 control_R2.fq -m Fastq2vcf -o outdir 
```
####Bam files to final result.
```Bash
python weseq.py -c config.cfg --inbam case.bam -m Bam2vcf -o outdir 
```
####Bam files to final result in somatic.
```Bash
python weseq_somatic.py -c config.cfg --inbam case.bam,control.bam -m Bam2vcf -o outdir 
```

##Mapper in iseq
###BWA
BWA:Fast, accurate, memory-efficient aligner for short and long sequencing reads
BWA (Burrows-Wheeler Aligner) is an aligner using the Burrows-Wheeler transform to index the reference genome, which decreases memory usage compared to aligners using k-mer hashing.
BWA includes two read alignment algorithms, the first is usually meant when simply the "BWA algorithm" is mentioned. It is callable via the command bwa align. The second algorithm is "BWA-SW", it can be called via the command bwa bwasw. That tool is described in its own article BWA-SW.
###Bowtie/Bowtie2
Bowtie:An ultrafast, memory-efficient short read aligner.
Bowtie aligns short DNA sequences (reads) to the human genome at a rate of 25 million reads per hour on a typical workstation with 2 gigabytes of memory. Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory footprint small: 1.3 GB for the human genome. It supports alignment policies equivalent to Maq and SOAP but is much faster: about 35x faster than Maq and over 350x faster than SOAP when aligning to the human genome.
###TopHat2
TopHat:A fast splice junction mapper for RNA-Seq reads.
TopHat aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner BOWTIE, and then analyzes the mapping results to identify splice junctions between exons.
###STAR
STAR:Ultrafast universal RNA-seq aligner
MOTIVATION: Accurate alignment of high-throughput RNA-seq data is a challenging and yet unsolved problem because of the non-contiguous transcript structure, relatively short read lengths and constantly increasing throughput of the sequencing technologies. Currently available RNA-seq aligners suffer from high mapping error rates, low mapping speed, read length limitation and mapping biases.
##Variant Caller in iseq
###HaplotypeCaller
Call germline SNPs and indels via local re-assembly of haplotypes
###UnifiedGenotyper
Call SNPs and indels on a per-locus basis
###MuTect2
Call somatic SNPs and indels via local re-assembly of haplotypes
###Varscan2
VarScan is a tool that detects variants (SNPs and indels) in next-generation sequencing data. The new release (VarScan 2) is implemented in Java, and includes several new features.
###LoFreq
LoFreq-Star (also known as LoFreq or LoFreq version 2) is a fast and sensitive variant-caller for inferring single-nucleotide variants (SNVs) from high-throughput sequencing data. It is designed to robustly call low-frequency variants in next-gen sequencing data-sets. LoFreq has been used to call rare variants in viral and bacterial sequencing datasets and can be used to study mitochondrial heteroplasmy and rare somatic mutations in heterogeneous tumors.
##Authors
The authors list is maintained in the AUTHORS file. See also the Contributors list at github.

##License
Licensed under the BSD License. See the [LICENSE.txt](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT) file.
