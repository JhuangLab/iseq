# iseq:A convenient integreted analysis tool of NGS Variant Discovery
iseq is an integrated analysis pipeline for NGS panel sequencing data. If you have any question about this tool, please contact us
(lee_jianfeng@sjtu.edu.cn)

## Feature

- Multiple sequence alignment softwares
- Multiple variantion detection softwares
- Easy to use
- Extend easily
- Single configuration file do all things
- Output log be divided into in individual file

## Requirements
- Java
- R
- Git
- GATK
- TVC
- Bwa
- STAR
- Bowtie
- Bowtie2
- Tophat2
- samtools
- Picard
- ANNOVAR
- annovarR
- BioInstaller


## Install
### pip
```bash
pip install iseq
```

### source code
```bash
# Download the source code and gunzip
python setup.py install
# Or use pip
pip install .
```

## Usage

**Without paired normal control sample**

```bash
# fastq2vcf mode
panel.py -c config.cfg \
         -s A01A \
         -m fastq2vcf \
         -1 A01A_1.fq.gz \
         -2 A01A_2.fq.gz \
         --bamprocess 00101111 \
         -o outdir

# fastq2bam
panel.py -c config.cfg \
         -s A01A \
         -m fastq2bam \
         -1 A01A_1.fq.gz \
         -2 A01A_2.fq.gz \
         --bamprocess 00101111 \
         -o outdir

# bam2vcf
panel.py -c config.cfg \
         -s A01A \
         -m bam2vcf \
         --in_bam A01A.bam \
         --bamprocess 00000000 \
         -o outdir

# genomeindex mode
panel.py -c config.cfg -m genomeindex
```

**With paired normal control sample**

```bash
# fastq2vcf mode
panel_somatic.py -c config.cfg \
                 -s A01 \
                 -m fastq2vcf \
                 -1 A01A_1.fq.gz \
                 -2 A01A_2.fq.gz \
                 -3 A01C_1.fq.gz \
                 -4 A01C_2.fq.gz \
                 --bamprocess 00101111 \
                 -o outdir

# fastq2bam mode
panel_somatic.py -c config.cfg \
                 -s A01 \
                 -m fastq2bam \
                 -1 A01A_1.fq.gz \
                 -2 A01A_2.fq.gz \
                 -3 A01C_1.fq.gz \
                 -4 A01C_2.fq.gz \
                 --bamprocess 00101111 \
                 -o outdir

# bam2vcf mode
panel_somatic.py -c config.cfg \
                 -s A01 \
                 -m bam2vcf \
                 --case_in_bam A01A.bam \
                 --control_in_bam A01C.bam \
                 --bamprocess 00000000 \
                 -o outdir

# genomeindex mode
panel_somatic.py -c config.cfg -m genomeindex
```
