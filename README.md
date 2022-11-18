# te_hic

# Why? 

Analysis of TEs in HiC data. 

One of the first steps in HiC analysis is to discard multimapping reads, and low quality mapping reads.
We beleive this is a mistake when you consider that most of the genome consists of badly mapping 
repeats. 

This is a set of scripts that aims to fill that gap by looking at HiC data from the perspective of TEs.

See the Hutchins lab here:
https://www.chrom-lab.org/

## What? 

Work in progress.

Doubly so, I as I rework the entire entry platform

## How?

Installation:

stick bin on your PATH:

export PATH=/path/to/te_hic/bin:$PATH

Requires:
python3
bowtie2
samtools
pysam


To fix: 
You will also need to premake the genome annotation files. Instructions are given for this in 
te_hic/genomes/

You need to run the scripts in there, and make the .glb files which will be needed for later steps

Workflow:

te_hic consists of a series of wrapper scripts around python scripts that are submitted to a cluster.

Todo!

All should go through te_hic command now.

## PBS-script workflow 

If you are using a PBS-torque cluster environment

A typical run:
1. setup a 'samples' folder with a structure like this:
    samples/
    samples/sample_name1/
                         FASTQ_seq1_1.fq.gz # PE 1
                         FASTQ_seq1_2.fq.gz # PE 2
                         FASTQ_seq2_1.fq.gz # PE 1
                         FASTQ_seq2_2.fq.gz # PE 2
    samples/sample_name2/
                         FASTQ_seq1_1.fq.gz # PE 1
                         FASTQ_seq1_2.fq.gz # PE 2
    and so on...

As HiC data is so large it often comes split up. The directory information is important as it will tell 
te_hic when to merge the FASTQ data. Biological replicates should be at the 'sample_name' level. Whilst technical 
replicates (i.e. extra sequencing of the same sample, or split lanes from the sequencing run) should be within the 
sample_name folder (e.g. FASTQ_seq1 and FASTQ_seq2 in the 'sample_name1' exampe above).

te_hic will merge the technical replicates, but will keep the biological replicates. 

2. Step 1. Align the reads using bowtie2.  
te_hic.1.align <qstat -q queue name> <species>
valid species names are hg38, mm10 more could be made. See the /te_hic/genome/ folder for instructions.
The qstat q name should be a valid q on your PBS cluster to submit the jobs to

2. Step 2. collect valid FASTQ pairs:
te_hic.2.get_valid_pairs <qstat -q queue name>

3. Step 3. Annotate reads:

## SLURM-script workflow 

If you are using a SLURM cluster environment

Todo?

## License

Release License:
MIT license:
Copyright (C) 2019-2022 Andrew Hutchins
    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
Except as contained in this notice, the name(s) of the above copyright holders shall not be used in advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization.
    
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


