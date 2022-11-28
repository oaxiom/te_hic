# te_hic

Analysis of TEs in HiC data. 

## Why? 

One of the first steps in HiC analysis is to discard multimapping reads,  low quality mapping reads, or 
otherwise disregard repeat regions of the genome.

We beleive this is a mistake when you consider that most of the genome consists of badly mapping 
repeats. 

This is a set of scripts that aims to fill that gap by looking at HiC data from the perspective of TEs.

See the Hutchins lab here:
https://www.chrom-lab.org/

## What? 

te_hic is designed to assist in the analysis of TEs in HiC data. 

It differs from most analysis software as it considers multiple mapped reads and does not delete 
repeat regions from the genome. It is primarily targeted at the analysis of LINEs, SINEs, LTRs,
DNA and retroposons (in human). The analysis pipeline excludes simple repeats, low complexity and
satellites. (The last of these is potentially very interesting, but the genome annotations
are not great for these, so we omit).


## How?

Installation:

Stick te_hic/bin on your PATH:
export PATH=/path/to/te_hic/bin:$PATH

Requires:
python3 (>3.7)
pysam

### Step 0: Build the genome indices for te_hic

Build the genome annotations.

This is not currently automated, but to do you need to go into te_hic and execute these commands:

```
cd te_hic/genome
./download_data.sh 

# hg38 genome
python make_hg38.py
python genome_freq_hg38.py

# mm10 genome:
TODO!

```

### Step 1: Align to the genome

Align your FASTQ HiC data to the genome using your aligner of choice.

We suggest something like :

```
bowtie2 --end-to-end --very-sensitive -x species_index -U read1.fq.gz | samtools view -b | samtools sort -n >${out}.p1.bam
bowtie2 --end-to-end --very-sensitive -x species_index -U read2.fq.gz | samtools view -b | samtools sort -n >${out}.p2.bam
```

You may need to merge the resulting BAM files. At the end there should be a single bam file
per read pair per sample.

### Step 2: Use te_hic

The 

## License

Release License:
MIT license:
Copyright (C) 2019-2022 Andrew Hutchins
    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
Except as contained in this notice, the name(s) of the above copyright holders shall not be used in advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization.
    
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


