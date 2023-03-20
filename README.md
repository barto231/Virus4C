# 4Cvirus

## Prepare References

You need 3 references to run the viral integration portion of the code. You will need the 3 following references:

1. The host reference genome (Ex. hg19, hg38, mm9)
2. The viral reference
3. A concatenated file of the host reference genome and the virusal reference genome

Each reference needs to be indexed and with bwa and samtools using the following:
```
bwa index host_reference.fa

samtools faidx host_reference.fa
```
```
bwa index viral_reference.fa

samtools faidx viral_reference.fa
```
```
cat host_reference.fa viral_reference.fa > host_reference+viral_reference.fa

bwa index host_reference+viral_reference.fa

samtools faidx host_reference+viral_reference.fa
```

## Compiling


The required external libraries (downloaded with the source code) must be compiled with
```
./build_libs.sh
```

Then, run
```
cmake -DCMAKE_BUILD_TYPE=Release . && make
```
## Required Software
- A Unix like shell (e.g. Bash v3.2+)
- Python
- Python Libraries
  - NumPy
  - PyFaidx
  - PySam
- Bowtie2 v2.3+ =
- SAMtools v1.3+ 
- sdust (https://github.com/lh3/sdust)
- R v3.5+ 
- The following R packages available from CRAN:
  - optparse
  - caTools
  - config
- The following R packages available from Bioconductor:
  - shortRead
  - genomicRanges
  - genomicAlignments
  - BSgenome
  - BSgenome of interest
  - rtracklayer (only if you want to generate bigWig files)
- The peakC package available from https://github.com/deWitLab/peakC/.

## Files required to run the pipepline


