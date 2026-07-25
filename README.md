# CSHA 2026 : Advanced Technologies for 3D Genome Research
## Suzhou, China / July 19 - 26, 2026
**Minji Kim**, University of Michigan

**Yijun Ruan**, Zhejiang University

**Fuchou Tang**, Peking University

This program offers a deep dive into the latest advancements in single-cell 3D genomics, tailored for graduate students, postdoctoral fellows, and early-career scientists worldwide. It is designed for participants with training in biology, molecular biology, bioinformatics, genetics, and related disciplines. The curriculum centers on tutorial-style instruction covering the principles and workflows of single-cell 3D genome profiling, complemented by oral teaching and curated video clips that walk through key experimental steps and practical considerations.

In parallel, the program features a dry-lab workshop that guides participants through core analytical skills, including data processing, quality control, visualization, and integrative interpretation of single-cell 3D genomics datasets. Expert-led lectures and discussions are integrated throughout the week to provide a coherent conceptual framework and up-to-date perspectives in the field.
##
## Hands-on Workshop: Computational Analysis of 3D Genomic Data
## Installation & Quick Start
```bash
micromamba create -n CSHA
micromamba activate CSHA

micromamba install \
-c conda-forge -c bioconda -c r \
--strict-channel-priority -y \
matplotlib biopython cooler sra-tools bwa minimap2 samtools pigz chromap "pairtools==0.3" "numpy<=1.23"
pip install runHiC

micromamba install \
-c conda-forge -c bioconda -c r \
--strict-channel-priority -y \
cooltools ucsc-bedgraphtobigwig pomegranate

pip install TADLib RaichuNorm statsmodels
pip install -U hicpeaks
```
## Processing Hi-C Data with runHiC
<img width="798" height="364" alt="image" src="https://github.com/user-attachments/assets/624c017d-85b2-4967-a05f-7f9e54c7286b" />

### Data Preparation
```bash
mkdir ref
mkdir ref/hg38

mkdir fastq
mkdir fastq/HiC-SRA
mkdir fastq/HiC-FASTQ
mkdir fastq/HiC-gzip

mkdir result
mkdir result/runHiC
mkdir result/runHiC/CSHA_20260723

mkdir command

cd fastq/HiC-SRA

prefetch -o SRR027956.sra SRR027956
prefetch -o SRR027958.sra SRR027958

for i in ./*.sra; do fastq-dump --split-files $i; done
for i in ./*.fastq; do gzip -c $i > `basename $i`.gz; done

cd ..
mv ./HiC-SRA/*.fastq ./HiC-FASTQ
mv ./HiC-SRA/*.gz ./HiC-gzip

cd ../ref
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

cd hg38
cp /data3/psg/PIPEref/hg38/* .

cd fastq
vi datasets.tsv

SRR027956 GM06990 R1 HindIII
SRR027958 GM06990 R2 HindIII
```
```bash
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes

awk 'BEGIN{OFS="\t"}
$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
    print $1, $2
}' ../ref/hg38.chrom.sizes \
> ../ref/hg38.autosomes.chrom.sizes

sort -k1,1V ../ref/hg38.autosomes.chrom.sizes -o ../ref/hg38.autosomes.chrom.sizes

awk 'BEGIN{OFS="\t"}
NF >= 2 && $2 ~ /^[0-9]+$/ {
    print $1, 0, $2
}' ../ref/hg38.autosomes.chrom.sizes \
> ../ref/hg38.autosomes.view.bed

sort -k1,1V ../ref/hg38.autosomes.view.bed -o ../ref/hg38.autosomes.view.bed

awk 'BEGIN{OFS="\t"}
$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
    print $1, $2
}' ../ref/mm10.chrom.sizes \
> ../ref/mm10.autosomes.chrom.sizes

sort -k1,1V ../ref/mm10.autosomes.chrom.sizes -o ../ref/mm10.autosomes.chrom.sizes

awk 'BEGIN{OFS="\t"}
NF >= 2 && $2 ~ /^[0-9]+$/ {
    print $1, 0, $2
}' ../ref/mm10.autosomes.chrom.sizes \
> ../ref/mm10.autosomes.view.bed

sort -k1,1V ../ref/mm10.autosomes.view.bed -o ../ref/mm10.autosomes.view.bed
```
### Mapping
```bash
runHiC mapping \
-m ../fastq/datasets.tsv \
-p ../ref/ \
-g hg38 \
-f ../fastq/HiC-gzip \
-F FASTQ \
-A bwa-mem \
-i ../ref/hg38/hg38.fa \
-t 20 \
--include-readid \
--drop-seq \
--logFile runHiC-mapping.log
```
### Filtering
```bash
runHiC filtering \
-m ../fastq/datasets.tsv \
--pairFolder pairs-hg38/ \
--logFile runHiC-filtering.log \
--nproc 20
```
### Binning
```bash
runHiC binning \
-f filtered-hg38/ \
--logFile runHiC-binning.log \
--nproc 20
```
### Moving
```bash
mv *hg38* ../result/runHiC/CSHA_20260723/.
mv *log ../result/runHiC/CSHA_20260723/.
```
### Results
```bash
GM06990-HindIII-allReps-filtered.mcool	101M
GM06990-HindIII-R1-filtered.mcool		56M
GM06990-HindIII-R2-filtered.mcool		65M
```
