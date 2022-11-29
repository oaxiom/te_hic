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

te_hic is all packed in a single entry point now. To use, supply the aligned read1 and read2 bams
and the genome. Other options are available to control the q read threshold, minimum distance to consider
and resolutions of arrays to build.

```
te_hic -1 ${out}.p1.bam -2 ${out}.p2.bam -g hg38 -l sample_label
```

These are the full optinos for te_hic:

```

usage: te_hic [-h] [-l LABEL] [-q QUAL] [-d MINDIST] [-r RESOLUTIONS [RESOLUTIONS ...]] -g GENOME -1 READ1 -2 READ2

HiC data analysis, preserving TE information

required arguments:
  -g GENOME, --genome GENOME
                        Genome assembly to use, valid genomes: {'hg38'}
  -1 READ1, --read1 READ1
                        the BAM alignment file containing the aligned reads pair 1
  -2 READ2, --read2 READ2
                        the BAM alignment file containing the aligned reads pair 2

optional arguments:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        label for the sample name and output files, default=te_hic
  -q QUAL, --qual QUAL  q threshold for read quality filtering, default=10
  -d MINDIST, --mindist MINDIST
                        Minimum distance (in base pairs), default=5000
  -r RESOLUTIONS [RESOLUTIONS ...], --resolutions RESOLUTIONS [RESOLUTIONS ...]
                        Default resolutions for the matrices default=300 150 50

Minimal usage usage: te_hic -a read1.bam -b read2.bam -g genome

```

## HiC example

A toy example is included in te_hic/test/HiC-analysis (the 'HiC-analysis % ' is the user prompt, type
in the commands after the % to get it it run)

```
HiC-analysis % te_hic -1 SRR1030718.2_1.100k.bam -2 SRR1030718.2_2.100k.bam -g hg38 -l toy_example

INFO    : Arguments:
INFO    :   Read1: SRR1030718.2_1.100k.bam
INFO    :   Read2: SRR1030718.2_2.100k.bam
INFO    :   Genome: hg38
INFO    :   Quality thresold: 10 (default is 10)
INFO    :   Minimum contact distance: 5000 (default is 5000)
INFO    :   Matrix resolutions to build: [300, 150, 50] kbp (default is [300, 150, 50])
INFO    : Loaded '/Users/andrew/Tools/te_hic/tehiclib/../genome/hg38_glb_gencode_tes.glb' binary file with 5665930 items
INFO    : Loaded '/Users/andrew/Tools/te_hic/tehiclib/../genome/hg38_te_genome_freqs.glb' binary file with 1058 items
INFO    : Stage 1: Collect valid read pairs
INFO    : Started SRR1030718.2_1.100k.bam and SRR1030718.2_2.100k.bam
[E::idx_find_and_load] Could not retrieve index file for 'SRR1030718.2_1.100k.bam'
[E::idx_find_and_load] Could not retrieve index file for 'SRR1030718.2_2.100k.bam'
INFO    : 
collect_valid_pairs() stats:
INFO    :   Aligned:
INFO    :     Reads processed           : 99,970
INFO    :     Correctly paired          : 85,809 (85.83%)
INFO    :   Rejected reads:
INFO    :     Low quality               : 4,213 (4.21%)
INFO    :     One pair aligned          : 12,252 (12.26%)
INFO    :     Not canonical chromosome  : 312 (0.31%)
INFO    :     No pairs aligned          : 1,909 (1.91%)
INFO    :     Duplicates                : 25 (0.05%)
INFO    :   Rejected reads (by criteria):
INFO    :     Too close                 : 28,649 (28.66%)
INFO    :   Final:
INFO    :     Kept reads                : 52,610 (52.63%)
INFO    :     Kept short-range (<20kb)  : 2,421 (2.42%)
INFO    :     Kept long-range (>20kb)   : 50,214 (50.23%)
INFO    : Intermediate file: Saved 52,610 pairs
INFO    : Stage 2: Assign to hg38 genome features
Processed 52,610 reads
INFO    : Intermediate file: Pair assignments 52,610 saved
INFO    : Stage 3: Quantify links
INFO    : Measures anchors...
INFO    :   TE <-> TE : 15,578 (29.61%)
INFO    :   TE <-> -- : 25,297 (48.08%)
INFO    :   -- <-> -- : 11,735 (22.31%)
INFO    : Saved 'stage3.te_hic_te-nn_anchor_frequencies.tsv'
INFO    : Stage 4: Build Matrices
INFO    : Building in-memory matrices for resolution 300000 kbp
INFO    : Saved BED bins: "matrices_te_hic/300000/te_hic_300000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/300000/te_hic_300000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/300000/te_hic_300000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/300000/te_hic_300000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/300000/te_hic_300000.nnnn.raw.matrix"
INFO    : Building in-memory matrices for resolution 150000 kbp
INFO    : Saved BED bins: "matrices_te_hic/150000/te_hic_150000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/150000/te_hic_150000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/150000/te_hic_150000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/150000/te_hic_150000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/150000/te_hic_150000.nnnn.raw.matrix"
INFO    : Building in-memory matrices for resolution 50000 kbp
INFO    : Saved BED bins: "matrices_te_hic/50000/te_hic_50000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/50000/te_hic_50000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/50000/te_hic_50000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/50000/te_hic_50000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/50000/te_hic_50000.nnnn.raw.matrix"

HiC-analysis % ls -l
total 27024
-rw-rw-r--  1 user  user  4208603 Nov 23 10:56 SRR1030718.2_1.100k.bam
-rw-rw-r--  1 user  user  4318275 Nov 23 11:03 SRR1030718.2_2.100k.bam
drwxr-xr-x  6 user  user      192 Nov 28 16:35 matrices_te_hic
-rw-r--r--  1 user  user   961678 Nov 29 07:32 stage1.int.te_hic.tsv.gz
-rw-r--r--  1 user  user  1111507 Nov 29 07:32 stage2.int.te_hic.tsv.gz
-rw-r--r--  1 user  user       93 Nov 29 07:32 stage3.te_hic_crude_measures.txt
-rw-r--r--  1 user  user    79178 Nov 29 07:32 stage3.te_hic_te-nn_anchor_frequencies.tsv
-rw-r--r--  1 user  user  1357255 Nov 29 07:32 stage3.te_hic_te-te_anchor_frequencies.tsv
-rw-r--r--  1 user  user       71 Nov 28 13:14 test_script.sh

HiC-analysis % ll matrices_te_hic
total 0
drwxr-xr-x  8 user  user  256 Nov 28 16:36 150000
drwxr-xr-x  8 user  user  256 Nov 28 16:36 300000
drwxr-xr-x  8 user  user  256 Nov 28 16:36 50000

HiC-analysis % ll matrices_te_hic/150000 
total 4760
-rw-r--r--  1 andrew  staff  585851 Nov 29 07:32 te_hic_150000.all.raw.matrix
-rw-r--r--  1 andrew  staff  130041 Nov 29 07:32 te_hic_150000.nnnn.raw.matrix
-rw-r--r--  1 andrew  staff  287127 Nov 29 07:32 te_hic_150000.tenn.raw.matrix
-rw-r--r--  1 andrew  staff  181403 Nov 29 07:32 te_hic_150000.tete.raw.matrix
-rw-r--r--  1 andrew  staff  603430 Nov 29 07:32 te_hic_150000_abs.bed

HiC-analysis % wc *.tsv
     904    5428   79178 stage3.te_hic_te-nn_anchor_frequencies.tsv
   11570   69423 1357255 stage3.te_hic_te-te_anchor_frequencies.tsv
   12474   74851 1436433 total
   
HiC-analysis % head stage3.te_hic_te-nn_anchor_frequencies.tsv
name	count	genome_count	genome_percent	RPM	RPM per kbp of TE
Alu:Alu:SINE	7	239798	0.007743788326674956	133.05455236647026	0.5548609761819125
AluJb:Alu:SINE	471	31276949	1.0100254070518015	8952.670594943927	0.2862386160153897
AluJo:Alu:SINE	291	17372139	0.5609978698636967	5531.267819806121	0.31839877747962536
AluJr4:Alu:SINE	79	4637982	0.1497741885709162	1501.6156624215928	0.32376487498692164
AluJr:Alu:SINE	301	19812811	0.6398144043754208	5721.345751758221	0.28877001611524084
AluSc5:Alu:SINE	28	1850062	0.05974398668556419	532.218209465881	0.2876758776007945
AluSc8:Alu:SINE	88	6125921	0.19782415003433296	1672.6858011784832	0.27305050149658855
AluSc:Alu:SINE	125	9718052	0.31382470927872713	2375.9741494012546	0.2444907836880534
AluSg4:Alu:SINE	19	1990659	0.06428428063032403	361.1480707089907	0.18142136383428337

```

## License

Release License:
MIT license:
Copyright (C) 2019-2022 Andrew Hutchins
    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
Except as contained in this notice, the name(s) of the above copyright holders shall not be used in advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization.
    
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


