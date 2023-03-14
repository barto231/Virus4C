# 4Cvirus

# Prepare References

You need 3 references to run the viral integration portion of the code. You will need the 3 following references:

1. The host reference genome (Ex. hg19, hg38, mm9)
2. The viral reference
3. A concatenated file of the host reference genome and the virusal reference genome

Each reference needs to be indexed and with bwa and samtools using the following:

bwa index host_reference.fa

samtools faidx host_reference.fa


bwa index viral_reference.fa

samtools faidx viral_reference.fa


cat host_reference.fa viral_reference.fa > host_reference+viral_reference.fa

bwa index host_reference+viral_reference.fa

samtools faidx host_reference+viral_reference.fa


