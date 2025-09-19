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
```
data/bam/
    |-- reads_father.bam
    |-- reads_mother.bam
    |-- reads_mother.bam.bai
    \-- reads_son.bam
```
## Exit the Samtools container
```
exit 
```
## Workflow: Call variants with GATK HaplotypeCaller

Pull the GATK container
```
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```
Spin up the GATK container interactively
```
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

 Run the variant calling command
We need to provide the BAM input file (-I) as well as the reference genome (-R), a name for the output file (-O) and a list of genomic intervals to analyze (-L).

However, we don't need to specify the path to the index file; the tool will automatically look for it in the same directory, based on the established naming and co-location convention. The same applies to the reference genome's accessory files (index and sequence dictionary files, *.fai and *.dict).
```
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

The output file reads_mother.vcf is created inside your working directory in the container, so you won't see it in the VS Code file explorer unless you change the output file path. However, it's a small test file, so you can cat it to open it and view the contents. If you scroll all the way up to the start of the file, you'll find a header composed of many lines of metadata, followed by a list of variant calls, one per line.


```
exit
```

Where does NEXTFLOW fit into all of these?

In order to make these stepwise commands, we can create a pipeline into a two-step workflow that uses containers to execute the work.


## Writing a single-stage workflow that runs Samtools index on a BAM file

```
genomics-1.nf
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

         // Save output files to the directory given by `params.outdir`
        // Instead of copying, create symbolic links
    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """
}
```

Even though the data files we're using here are very small, in genomics they can get very large. For the purposes of demonstration in the teaching environment, we're using the 'symlink' publishing mode to avoid unnecessary file copies. You shouldn't do this in your final workflows, since you'll lose results when you clean up your work directory.

Then we add the workflow at the end since we need to set up a channel to feed the input to the SAMTOOLS_INDEX process; then we can call the process itself to run on the contents of that channel.

```
workflow {

    // Create input channel (single file via CLI parameter)
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
}

```
Then run using:

```
nextflow run genomics-1.nf
```






