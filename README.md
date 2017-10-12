# iseq

iseq is a integrated analysis toolkits and pipeline for NGS panel sequencing data. If you have any question about this tool, please contact [us](lee_jianfeng@sjtu.edu.cn).

**Python class of iseq:** ReffaFile, FastqFile, BamFile, SamFile, VcfFile, CsvFile, MpileupFile and ResultFile.

**Processor of iseq:** preprocess, variantcaller, refinement

**Pipeline:** panel, panel_somatic

## Feature

- Multiple sequence alignment softwares
- Multiple variantion detection softwares
- Easy to use
- Extend easily
- Single configuration file do all things
- Output log be divided into in individual file

## Install

### source code

```bash
# Download the source code and gunzip
./install.R
```

## Usage

[configuration file](https://github.com/JhuangLab/iseq/blob/master/data/config.cfg) is an important file to run iseq. Three type of parameters can be found in this file:

- tools
- extra files
- non-file type parameters

**Tools:** Command line tools, class of iseq can use these tools by `self.cfg`, e.g `self.cfg["gatk"]` (Need to be download and install)

**Extra files:** Some of files required to run tools (Need to be download)

**Non-file type parameters:** Other parameters (Do not need to be download anything)

### Germline Mode

```bash
# fastq2vcf mode
panel -c config.cfg \
         -s A01A \
         -m fastq2vcf \
         -1 A01A_1.fq.gz \
         -2 A01A_2.fq.gz \
         --bamprocess 00101111 \
         -o outdir

# fastq2bam
panel -c config.cfg \
         -s A01A \
         -m fastq2bam \
         -1 A01A_1.fq.gz \
         -2 A01A_2.fq.gz \
         --bamprocess 00101111 \
         -o outdir

# bam2vcf
panel -c config.cfg \
         -s A01A \
         -m bam2vcf \
         --in_bam A01A.bam \
         --bamprocess 00000000 \
         -o outdir

# genomeindex mode
panel -c config.cfg -m genomeindex
```

### Somatic Mode

```bash
# fastq2vcf mode
panel_somatic -c config.cfg \
                 -s A01 \
                 -m fastq2vcf \
                 -1 A01A_1.fq.gz \
                 -2 A01A_2.fq.gz \
                 -3 A01C_1.fq.gz \
                 -4 A01C_2.fq.gz \
                 --bamprocess 00101111 \
                 -o outdir

# fastq2bam mode
panel_somatic -c config.cfg \
                 -s A01 \
                 -m fastq2bam \
                 -1 A01A_1.fq.gz \
                 -2 A01A_2.fq.gz \
                 -3 A01C_1.fq.gz \
                 -4 A01C_2.fq.gz \
                 --bamprocess 00101111 \
                 -o outdir

# bam2vcf mode
panel_somatic -c config.cfg \
                 -s A01 \
                 -m bam2vcf \
                 --case_in_bam A01A.bam \
                 --control_in_bam A01C.bam \
                 --bamprocess 00000000 \
                 -o outdir

# genomeindex mode
panel_somatic -c config.cfg -m genomeindex
```
