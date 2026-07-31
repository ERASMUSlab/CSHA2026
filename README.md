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

## Identifying 3D chromatin structure
### assessing the reproducibility of Hi-C data with HiC-Rep
<img width="875" height="457" alt="image" src="https://github.com/user-attachments/assets/538821ae-fd93-4b65-86f7-496b326883c8" />

### Data Preparation
```bash
mkdir input_data
mkdir input_data/hicrep
mkdir result/hicrep

cd input_data/hicrep

wget https://github.com/ERASMUSlab/CSHA2026/releases/download/ver1/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/ver1/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/ver1/K562_insituhic_rep1_4DNFI9G9FRJJ.ice.mcool
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/ver1/K562_insituhic_rep2_4DNFIRHH2E7D.ice.mcool
```
### CLI Analysis
```bash
cd ../../command/

hicrep \
../input_data/hicrep/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool \
../input_data/hicrep/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool \
../result/hicrep/GM12878_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM

hicrep \
../input_data/hicrep/K562_insituhic_rep1_4DNFI9G9FRJJ.ice.mcool \
../input_data/hicrep/K562_insituhic_rep2_4DNFIRHH2E7D.ice.mcool \
../result/hicrep/K562_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM

hicrep \
../input_data/hicrep/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool \
../input_data/hicrep/K562_insituhic_rep1_4DNFI9G9FRJJ.ice.mcool \
../result/hicrep/GM12878rep1_K562rep1_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM

hicrep \
../input_data/hicrep/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool \
../input_data/hicrep/K562_insituhic_rep2_4DNFIRHH2E7D.ice.mcool \
../result/hicrep/GM12878rep1_K562rep2_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM

hicrep \
../input_data/hicrep/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool \
../input_data/hicrep/K562_insituhic_rep1_4DNFI9G9FRJJ.ice.mcool \
../result/hicrep/GM12878rep2_K562rep1_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM

hicrep \
../input_data/hicrep/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool \
../input_data/hicrep/K562_insituhic_rep2_4DNFIRHH2E7D.ice.mcool \
../result/hicrep/GM12878rep2_K562rep2_hicrep.txt \
--binSize 50000 \
--h 4 \
--dBPMax 500000 \
--excludeChr chrM
```
### Result
```bash
GM12878_hicrep.txt                 1023
GM12878rep1_K562rep1_hicrep.txt    1.1K
GM12878rep1_K562rep2_hicrep.txt    1.1K
GM12878rep2_K562rep1_hicrep.txt    1.1K
GM12878rep2_K562rep2_hicrep.txt    1.1K
K562_hicrep.txt                    1014
```

### API Analysis
```python
from hicrep import hicrepSCC
from hicrep.utils import readMcool
import numpy as np
import matplotlib.pyplot as plt

h = 4
dBPMax = 500000 
bDownSample = False
resolution = 50_000

f_hic1 = f"/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool"
f_hic2 = f"/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool"

cool1, binSize1 = readMcool(f_hic1, resolution)
cool2, binSize2 = readMcool(f_hic2, resolution)

scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, excludeChr="chrM")
sccSub = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, np.array(['chr16'], dtype=str))

f_hic_many = [
    "/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/GM12878_insituhic_rep1_4DNFIUOVQH68.ice.mcool",
    "/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/GM12878_insituhic_rep2_4DNFIIT7LQ6M.ice.mcool",
    "/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/K562_insituhic_rep1_4DNFI9G9FRJJ.ice.mcool",
    "/data3/psg/NGS_2026/CSHA_HiC/input_data/hicrep/K562_insituhic_rep2_4DNFIRHH2E7D.ice.mcool",
]

names = [
'GM12878 rep1',
'GM12878 rep2',
'K562 rep1',
'K562 rep2',
]

cools = [readMcool(f, resolution)[0] for f in f_hic_many]

from itertools import combinations

n = len(f_hic_many)
mat = np.eye(n)

for i, j in combinations(range(n), 2):

    scc = hicrepSCC(cools[i], cools[j], h, dBPMax, bDownSample,
                    np.array(['chr16'], dtype=str))
    
    mat[i, j] = mat[j, i] = scc[0]
```
### Result
```python
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='viridis')

ax.set_xticks(range(len(names)))
ax.set_yticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha='right')
ax.set_yticklabels(names)

thresh = (mat.max() + mat.min()) / 2
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, '{:.3g}'.format(mat[i, j]), ha='center', va='center',
                color='white' if mat[i, j] < thresh else 'black')

fig.colorbar(im, ax=ax, label='HiCRep (chr16)')
plt.tight_layout()
plt.show()
```
<img width="588" height="470" alt="image" src="https://github.com/user-attachments/assets/a65a572c-4fd4-431d-95d5-3ee70534b955" />


### A/B Compartments & TADs with cooltools
<img width="2091" height="845" alt="image" src="https://github.com/user-attachments/assets/cfdde319-d2a9-4442-9ab7-bf75dead04c3" />

### A/B Compartments Analysis
```bash
mkdir ../result/ABcompartment
mkdir ../result/ABcompartment/CSHA_20260723

cooltools genome binnify \
../ref/hg38.autosomes.chrom.sizes 50000 \
> ../result/ABcompartment/CSHA_20260723/ModifiedChromSize.txt

cooltools genome genecov \
../result/ABcompartment/CSHA_20260723/ModifiedChromSize.txt \
hg38 \
> ../result/ABcompartment/CSHA_20260723/ModifiedGeneDensity.bedGraph

cooltools eigs-cis \
--view ../ref/hg38.autosomes.view.bed \
--phasing-track ../result/ABcompartment/CSHA_20260723/ModifiedGeneDensity.bedGraph \
--out-prefix ../result/ABcompartment/CSHA_20260723/ABcompartment \
--bigwig \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered.mcool::resolutions/50000
```
### A/B Compartments Results
```bash
ABcompartment.cis.bw	    849K
ABcompartment.cis.lam.txt	2.0K
ABcompartment.cis.vecs.tsv	5.3M
```
### TADs Analysis
```bash
mkdir ../result/TADs
mkdir ../result/TADs/CSHA_20260723

cooltools insulation \
--view ../ref/hg38.autosomes.view.bed \
-p 20 \
-o ../result/TADs/CSHA_20260723/25000_InsulationScore.tsv \
--bigwig \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered.mcool::resolutions/25000 \
25000
```
### TADs Results
```bash
25000_InsulationScore.tsv             7.2M
25000_InsulationScore.tsv.25000.bw    566K
```

### Loops with HICCUPS/raichu
<img width="915" height="357" alt="image" src="https://github.com/user-attachments/assets/05782d78-0b03-4399-817f-5d8f276946dc" />

### Loops Analysis with HICCUPS
```bash
mkdir ../result/Loops
mkdir ../result/Loops/CSHA_20260723

cooler dump --table chroms \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered.mcool::resolutions/10000 \
> ../result/Loops/all.chrom.sizes

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ {print $1, $2}' \
../result/Loops/all.chrom.sizes \
| sort -k1,1V \
> ../result/Loops/autosome.chrom.sizes

cooler dump \
--table pixels \
--join \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered.mcool::resolutions/10000 \
| awk 'BEGIN{OFS="\t"}
NR==FNR { keep[$1]=1; next }
($1 in keep) && ($4 in keep) { print }' \
../result/Loops/autosome.chrom.sizes - \
| cooler load -f bg2 --assembly hg38 \
../result/Loops/autosome.chrom.sizes:10000 \
- ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool

cooler dump \
--table pixels \
--join \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R1-filtered.mcool::resolutions/10000 \
| awk 'BEGIN{OFS="\t"}
NR==FNR { keep[$1]=1; next }
($1 in keep) && ($4 in keep) { print }' \
../result/Loops/autosome.chrom.sizes - \
| cooler load -f bg2 --assembly hg38 \
../result/Loops/autosome.chrom.sizes:10000 \
- ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R1-filtered_10000.autosome.cool

cooler dump \
--table pixels \
--join \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R2-filtered.mcool::resolutions/10000 \
| awk 'BEGIN{OFS="\t"}
NR==FNR { keep[$1]=1; next }
($1 in keep) && ($4 in keep) { print }' \
../result/Loops/autosome.chrom.sizes - \
| cooler load -f bg2 --assembly hg38 \
../result/Loops/autosome.chrom.sizes:10000 \
- ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R2-filtered_10000.autosome.cool

cooler balance \
-p 20 --force \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool

cooler balance \
-p 20 --force \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R1-filtered_10000.autosome.cool

cooler balance \
-p 20 --force \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R2-filtered_10000.autosome.cool

pyHICCUPS \
-O ../result/Loops/CSHA_20260723/Result_ICE_Loops.bedpe \
-p ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool \
--clr-weight-name weight \
--pw 2 3 4 \
--ww 5 6 7 \
--maxww 10 \
--min-local-reads 25 \
--only-anchors \
--min-marginal-peaks 3 \
--maxapart 3000000 \
--nproc 20 \
--logFile hiccups_ICE.log
```
### Loops Analysis with raichu
```bash
raichu \
--cool-uri ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool \
-n raichu_weight \
--force \
--logFile raichu.log

pyHICCUPS \
-O ../result/Loops/CSHA_20260723/Result_Raichu_Loops.bedpe \
-p ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool \
--clr-weight-name raichu_weight \
--pw 2 3 4 \
--ww 5 6 7 \
--maxww 10 \
--min-local-reads 25 \
--only-anchors \
--min-marginal-peaks 3 \
--maxapart 3000000 \
--nproc 20 \
--logFile hiccups_raichu.log
```
### Compare loops
```bash
peak-plot \
-O ../result/Loops/CSHA_20260723/ICE_chr21_32000000_34000000.png \
-p ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool \
-I ../result/Loops/CSHA_20260723/Result_ICE_Loops.bedpe \
-C chr21 \
-S 32000000 \
-E 34000000 \
--clr-weight-name weight

peak-plot \
-O ../result/Loops/CSHA_20260723/Raichu_chr21_32000000_34000000.png \
-p ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-allReps-filtered_10000.autosome.cool \
-I ../result/Loops/CSHA_20260723/Result_Raichu_Loops.bedpe \
-C chr21 \
-S 32000000 \
-E 34000000 \
--clr-weight-name raichu_weight
```
### Result
<img width="960" height="388" alt="image" src="https://github.com/user-attachments/assets/0c078b41-8f89-4bca-9509-0c440b033f8e" />

### stripe with stripenn
<img width="1160" height="453" alt="image" src="https://github.com/user-attachments/assets/ea17b0fe-efbc-408c-8870-6167bf10c23c" />

### Data Preparation
```bash
cd ../input_data
mkdir stripenn
cd stripenn
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/stripenn_input/GM12878_insituhic_4DNFIXP4QG5B_5000.cool.part_aa
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/stripenn_input/GM12878_insituhic_4DNFIXP4QG5B_5000.cool.part_ab
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/stripenn_input/GM12878_insituhic_4DNFIXP4QG5B_5000.cool.part_ac
wget https://github.com/ERASMUSlab/CSHA2026/releases/download/stripenn_input/GM12878_insituhic_4DNFIXP4QG5B_5000.cool.sha256

cat GM12878_insituhic_4DNFIXP4QG5B_5000.cool.part_* > GM12878_insituhic_4DNFIXP4QG5B_5000.cool
shasum -a 256 -c GM12878_insituhic_4DNFIXP4QG5B_5000.cool.sha256
cooler info GM12878_insituhic_4DNFIXP4QG5B_5000.cool

pip install stripenn
```

### Analysis
```bash
cd ../../command/
mkdir ../result/stripe

stripenn compute \
--cool ../input_data/stripenn/GM12878_insituhic_4DNFIXP4QG5B_5000.cool \
--out ../result/stripe/GM12878 \
-k chr16

peak-plot \
-O ../result/stripe/GM12878/ICE_GM12878_chr16_13900000_14700000.png \
-p ../input_data/stripenn/GM12878_insituhic_4DNFIXP4QG5B_5000.cool \
-C chr16 \
-S 13900000 \
-E 14700000 \
--clr-weight-name weight

peak-plot \
-O ../result/stripe/GM12878/ICE_GM12878_chr16_46700000_48300000.png \
-p ../input_data/stripenn/GM12878_insituhic_4DNFIXP4QG5B_5000.cool \
-C chr16 \
-S 46700000 \
-E 48300000 \
--clr-weight-name weight
```

### Result
```bash
ICE_GM12878_chr16_13900000_14700000.png    99K
ICE_GM12878_chr16_46700000_48300000.png    246K
result_filtered.tsv                        8.6K
result_unfiltered.tsv                      241K
stripenn.log                               235
```
<img width="960" height="376" alt="image" src="https://github.com/user-attachments/assets/3169123a-ec44-4531-84a9-ef19917f9439" />

### Enhancing Hi-C data with HiCFoundation
<img width="1019" height="394" alt="image" src="https://github.com/user-attachments/assets/487ce388-a21a-4f5f-97e1-d007b6fdf006" />

### Analysis
```bash
mkdir ../result/TADs
mkdir ../result/TADs/CSHA_20260723

bash hitad_DS_code.sh \
25000 \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R1-filtered_25000.autosome.cool \
../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R2-filtered_25000.autosome.cool

hitad \
-d hitad_dataset.tsv \
-p 20 \
-O ../result/TADs/CSHA_20260723/hitad_25000.bed \
--removeCache \
--logFile hitad.log

tad-plot \
-O ../result/TADs/CSHA_20260723/TADvis.png \
--dpi 500 \
-p ../result/runHiC/CSHA_20260723/coolers-hg38/GM06990-HindIII-R2-filtered_25000.autosome.cool \
-T ../result/TADs/CSHA_20260723/hitad_25000.bed \
-C 21 \
-S 30000000 \
-E 36000000
```
### Result
<img width="960" height="251" alt="image" src="https://github.com/user-attachments/assets/f7f246b3-0ee3-473e-a82f-b2682920ab0d" />


