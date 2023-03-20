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

- Reads in (compressed) FASTQ format.
- Configuration file (conf.yml)

  
  <BR>
  
  
| Name            | Description                                                                                                                                                                                   |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| fragFolder      | Path to the folder containing the fragment end libraries of the reference genomes                                                                                                             |
| normalizeFactor | Reads mapped to the 4C fragment end library are normalized to account for sequencing depth according to the normalizeFactor                                                                   |
| enzymes         | Enzyme names used in the viewpoint file and their corresponding recognition motifs                                                                                                            |
| genomes         | Genome names used in the viewpoint file plus corresponding BSgenome packages                                                                                                                  |
| bowtie2         | Path to corresponding bowtie2 index of reference genome. The reference genome assembly used to generate the index should match to the reference genome that was used to generate the BSgenome (So the UCSC reference genomes for most BSgenome packages). |
| maxY            | Maximal Y value in local 4C cis plot                                                                                                                                                          |
| plotView        | Number of bp to plot around viewpoint in local 4C cis plot                                                                                                                                    |
| xaxisUnit       | X-axis unit (Mb, Kb or bp)                                                                                                                                                                    |
| plotType        | Plots will either be saved as PDF or PNG                                                                                                                                                      |
| binSize         | Genome bin size used in the genome plot                                                                                                                                                       |
| qualityCutoff   | Q-score. Trim 3′-end of all sequences using a sliding window as soon as 2 out of 5 nucleotides have quality encoding less than the Q-score. 0 = no trimming                                   |
| trimLength      | Trim reads to defined capture length from 3′-end. 0 = no trimming                                                                                                                             |
| minAmountReads  | Minimum required amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed                                                    |
| maxAmountReads  | Max amount of reads that will be used for the analysis. For 4C experiments only 1 million reads are required, more than 1 million reads will only increase the complexity if you also increase the amount of template used for the PCR. However when sequencing many reads per experiment the pipeline may crash at the trim step due to lack of memory. If this is the case reduce the maxAmountReads. In the future the pipeline will be updated to be memory efficieint in the trim step so just all reads can be used.|
| readsQuality    | Bowtie2 minimum required mapping quality score for mapped reads                                                                                                                               |
| mapUnique       | Extract uniquely mapped reads, based on the lack of XS tag                                                                                                                                    |
| cores           | Number of CPU cores for parallelization                                                                                                                                                       |
| wSize           | The running mean window size                                                                                                                                                                  |
| nTop            | Top fragment ends discarded for calculation of normalizeFactor                                                                                                                                |
| nonBlind        | Only keep non-blind fragments                                                                                                                                                                 |
| wig             | Create wig files for all samples                                                                                                                                                              |
| bigwig          | Create bigWig files for all samples                                                                                                                                                           |
| plot            | Create viewpoint coverage plot for all samples                                                                                                                                                |
| genomePlot      | Create genomeplot for all samples (only possible if analysis is “all” in vpFile)                                                                                                              |
| tsv             | Create tab separated value file for all samples                                                                                                                                               |
| bins            | Count reads for binned regions                                                                                                                                                                |
| mismatchMax     | The maximum number of mismatches allowed during demultiplexing                                                                                                                                |
| chr_random     | Do not include random chromosomes in the frag genome                                                                                                                                |
| chr_fix     | Do not include fix chromosomes in the frag genome  (Fix patches denoted by chr__fix represent changes to the existing sequence.)                                                                                                                              |
| chrUn     | Do not include Unknown chromosomes in the frag genome                                                                                                                                |
| chrM     | Do not include chromosome M in the frag genome                                                                                                                                |
| prefix     | The chromosome name prefix. For UCSC reference genomes use chr.                                                                                                                                 |
  
  **Table 1.** Description of parameters that need to be defined in the configuration file.
  
  <BR>
  
* Viewpoint file
  * Experiment specific parameters for each 4C-seq experiment are organized in a viewpoint file. Parameters in this file are stored in a tab-delimited format, with each row containing information for a separate experiment: 
  
| expname    | primer               | firstenzyme | secondenzyme | genome | vpchr | vppos    | analysis | fastq           |
|------------|----------------------|-------------|--------------|--------|-------|----------|----------|-----------------|
| mESC_Sox2  | GAGGGTAATTTTAGCCGATC | DpnII       | Csp6I        | mm9    | 3     | 34547661 | all      | index1.fastq.gz |
| mESC_Mccc1 | TTGCACCCGTCTTCTTGATC | DpnII       | Csp6I        | mm9    | 3     | 35873313 | cis      | index1.fastq.gz |

**Table 2.** Example of a viewpoint file in which two experiments are demultiplexed from the same FASTQ file based on their primer sequence.

<BR>
 

| Name              | Description                                                                                                                                                                                                                                                                                  |
|-------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| expname           | Unique experiment name                                                                                                                                                                                                                                                                       |
| primer            | Primer sequence (only include the sequence that will be included in the Illumina read, so exclude the Illumina adapter). Use the RE1 motif if the primer sequence untill the RE1 motif is already removed in the fastq file.                                                                                                                                                                       |
| firstenzyme       | First restriction enzyme name (nearest to reading primer)                                                                                                                                                                                                                                    |
| secondenzyme      | Second restriction enzyme name                                                                                                                                                                                                                                                               |
| genome            | Reference genome of interest                                                                                                                                                                                                                                                                 |
| vpchr             | The chromosome that contains the viewpoint. This information is only required to calculate the stats, if you have multiple integrations or if you do not know your VP location, just use 1.                                                                                                                                                                                                                                     |
| vppos             | Coordinate of viewpoint position. Any bp position within the VP can be used except the RE motifs. This information is only required to calculate the stats, if you have multiple integrations or if you do not know your VP location, just use 1.                                                                                                                                                                                |
| analysis          | The final output tables will contain all reads (all) or only the reads that have been mapped to the VP chromosome (cis). For most analysis cis is sufficient and the generated output files will be smaller and therefore easier to process on local computers                               |
| fastq             | Name of the FASTQ file                                                                                                                                                                                                                                                                       |
| spacer (optional) | Spacer length. Number of nt included as spacer in the primer to enable out of phase sequencing. Default = 0. The spacer sequence will not be used for demultiplexing. If the spacer sequence is used as a barcode include the sequence in the primer sequence and set the spacer length to 0 |

**Table 3.** Description of parameters that are required in the viewpoint file for processing a 4C-seq experiment.

| Name           | Description                                                                                                                                |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| vpFile*        | path to the viewpoint file                                                                                                                 |
| fqFolder*      | path to the folder containing the FASTQ files                                                                                              |
| outFolder*     | path to the output folder                                                                                                                  |
| confFile       | path to configuration file – default is conf.yml in folder containing the pipeline script                                                  |
| mismatchMax    | The maximum number of mismatches allowed during demultiplexing                                                                             |
| qualityCutoff  | Q-score. Trim 3′-end of all sequences using a sliding window as soon as 2 out of 5 nucleotides has quality encoding less than the Q-score  |
| trimLength     | Trim reads to defined capture length from 3′-end                                                                                           |
| minAmountReads | Minimum required amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed |
| readsQuality   | Bowtie2 minimum required mapping quality score for mapped reads                                                                            |
| mapUnique      | Extract uniquely mapped reads, based on the lack of XS tag                                                                                 |
| cores          | Number of CPU cores for parallelization                                                                                                    |
| wSize          | The running mean window size                                                                                                               |
| nTop           | Top fragment ends discarded for normalization                                                                                              |
| nonBlind       | Only keep non-blind fragments                                                                                                              |
| wig            | Create wig files for all samples                                                                                                           |
| bigwig         | Create bigWig files for all samples                                                                                                        |
| plot           | Create viewpoint coverage plot for all samples                                                                                             |
| genomePlot     | Create genomeplot for all samples (only possible if analysis is “all” in vpFile)                                                           |
| tsv            | Create tab separated value file for all samples                                                                                            |
| bins           | Count reads for binned regions                                                                                                             |

**Table 4.** Description of parameters that are recognized by the pipe4C.R script. * are required. 
