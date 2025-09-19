## Per-sample variant calling

First, I will be building a variant calling pipeline that applies GATK variant calling to individual sequencing samples. The Variant Calling pipeline aims to detect and identify variations in a genome sequence as compariosn to a reference genome, such as SNPs and indels.

#### Dataset
The dataset for this pipeline includes: A reference genome, 3 whole genome sequencing samples, a list of genomic intervals (annotation file). 

#### Workflow
1. Index a BAM input file with Samtools

# Pull the Samtools container
```
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

