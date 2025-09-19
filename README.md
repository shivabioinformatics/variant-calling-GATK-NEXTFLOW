# Per-sample variant calling

First, I will be building a variant calling pipeline that applies GATK variant calling to individual sequencing samples. The Variant Calling pipeline aims to detect and identify variations in a genome sequence as compared to a reference genome, such as SNPs and indels.

## Dataset
The dataset for this pipeline includes: A reference genome, 3 whole genome sequencing samples, and a list of genomic intervals (annotation file). 

## Workflow -- Set up 
1. Index a BAM input file with Samtools

### Pull the Samtools container
```
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

### Spin up the Samtools container interactively
```
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```
### Run the indexing command
```
samtools index /data/bam/reads_mother.bam
```
This process is fast, and immediately, BAM index files will be generated. 
data/bam/
├── reads_father.bam
├── reads_mother.bam
├── reads_mother.bam.bai
└── reads_son.bam



## Exit the Samtools container
```
exit 
```
