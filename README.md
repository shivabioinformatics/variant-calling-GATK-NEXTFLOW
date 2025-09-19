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

Nextflow:
Resource: https://sateeshperi.github.io/nextflow_varcal/nextflow/nextflow_channels

Channels are how Nextflow handles file management, allowing complex tasks to be split up, run in parallel, and reduces ‘admin’ required to get the right inputs to the right parts of the pipeline.
Channels connect processes via their inputs and outputs.
Channels can store multiple items, such as files (e.g., fastq files) or values.
The number of items a channel stores determines how many times a process will run using that channel as input. When the process runs using one item from the input channel, we will call that run a task.
Channels are asynchronous, which means that output data from a set of processes will not necessarily be output in the same order as they went in.
However, the first element into a queue is the first out of the queue (First in-First out). This allows processes to run as soon as they receive input from a channel. Channels only send data in one direction, from a producer (a process/operator), to a consumer (another process/operator).

**Nextflow distinguishes between two different kinds of channels:**
value channels
queue channels
Value channels
A value channel is bound to a single value.
A value channel can be used an unlimited number of times since its content is not consumed. This is also useful for processes that need to reuse input from a channel, for example, a reference genome sequence file that is required by multiple steps within a process, or by more than one process.
The value factory method is used to create a value channel. Values are put inside parentheses () to assign them to a channel.


Queue channel
Queue channels are a type of channel in which data is consumed (used up) to provide input for a process/operator.
Queue channels can be created in two ways:
As the outputs of a process.
A queue channel can be explicitly created using channel factory methods such as Channel.of or Channel.fromPath.
Queue channel factory
Channel factories are used to explicitly create channels. Channel factories are called using the Channel.<method> syntax and return a specific instance of a Channel.

Queue (consumable) channels can be created using the following channel factory methods:

Channel.of
Channel.fromList
Channel.fromPath
Channel.fromFilePairs
Channel.fromSRA


The **of** Channel factory
When you want to create a channel containing multiple values, you can use the channel factory Channel.of.
Channel.of allows the creation of a queue channel with the values specified as arguments, separated by a comma ,.
You can specify a range of numbers as a single argument using the Groovy range operator ... This creates each value in the range (including the start and end values) as a value in the channel. The Groovy range operator can also produce ranges of dates, letters, or time. More information on the range operator can be found here.

The **fromList** Channel factory
You can use the Channel.fromList method to create a queue channel from a list object.
Create a new file queue_list.nf; add the following and execute it using nextflow run queue_list.nf:

**Processes**
A process is the way Nextflow executes commands you would run on the command line or custom scripts.
Processes can be thought of as particular tasks or steps in a workflow, e.g. an alignment step (bwa mem) in variant-calling analysis.
Processes are independent of each other (don’t require another process to execute) and cannot communicate/write to each other.
It is the channels that pass the data from one process to another, and we do this by having the processes define input and output channels.

In Nextflow, the process definition starts with the keyword process, followed by the process name and finally the process body delimited by curly brackets {}.
The process body must contain a string which represents the command or, more generally, a script that is executed by it.
In order to run the process, we need to add it to a workflow block below the process.
The workflow scope starts with the keyword workflow, followed by an optional name and finally the workflow body delimited by curly brackets {}.

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






