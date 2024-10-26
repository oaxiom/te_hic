# te_hic

Analysis of TEs in HiC data. 

## Why? 

One of the first steps in HiC analysis is to discard multimapping reads, low quality mapping reads, or 
otherwise disregard repeat regions of the genome.

We believe this is a mistake when you consider that most of the genome consists of badly 
mapping repeats. 

This is a tool that aims to fill that gap by looking at HiC data from the perspective of TEs.

See the Hutchins lab here:
https://www.chrom-lab.org/

## What? 

te_hic is designed to assist in the analysis of TEs in HiC data. 

It differs from most HiC analysis software as it considers multiple mapped reads and does not delete 
repeat regions from the genome. It is primarily targeted at the analysis of LINEs, SINEs, LTRs,
DNA and retroposons (in human). The analysis pipeline excludes simple repeats, low complexity and
satellites. (The last of these is potentially very interesting, but the genome annotations
are not great, so we omit).

## How?

Installation:

Stick te_hic/bin on your PATH:

```
export PATH=/path/to/te_hic/bin:$PATH
```

Requires:
python3 (>3.7)
pysam

### Step 0: Build the genome indices for te_hic

Build the genome annotations.

```
% te_hic_genome -h
usage: te_hic_genome [-h] -g GENOME

Builds the indices for te_hic

options:
  -h, --help            show this help message and exit

required arguments:
  -g GENOME, --genome GENOME
                        Genome assembly to use, valid genomes: {'mm10', 'rn7', 'mm39', 'xenTro10', 'hg38', 'danRer11'}

Example usage: te_hic_genome -g hg38


```

For example:

```
te_hic_genome -g hg38
```

This assumes you have wget on your machine.

Genome currently supported:

Human: hg38 - GRCh38.p14 
Mouse: mm10 - NCBIm38
Mouse: mm39 - NCBIm39
Rat: rn7 - mRatBN7.2
Zebrafish: danRer11
Xenopus: xenTro10


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

te_hic is all packed in a single entry point. To use, supply the aligned read1 and read2 bams
and the genome. Other options are available to control the q read threshold, minimum distance to consider
and resolutions of arrays to build.

Note that te_hic uses one CPU core, and has a peak memory usage of about 11 Gb. 

```
te_hic -1 ${out}.p1.bam -2 ${out}.p2.bam -g hg38 -l sample_label
```

These are the full options for te_hic:

```

% te_hic -h
usage: te_hic [-h] [-l LABEL] [-q QUAL] [-d MINDIST] [-r RESOLUTIONS [RESOLUTIONS ...]] [--keep-intermediate-files] -g GENOME -1 READ1 -2 READ2

HiC data analysis, preserving TE information

required arguments:
  -g GENOME, --genome GENOME
                        Genome assembly to use, valid genomes: {'mm10', 'danRer11', 'rn7', 'hg38', 'xenTro10', 'mm39'}
  -1 READ1, --read1 READ1
                        the BAM alignment file containing the aligned reads pair 1
  -2 READ2, --read2 READ2
                        the BAM alignment file containing the aligned reads pair 2

options:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        label for the sample name and output files, default=te_hic
  -q QUAL, --qual QUAL  q threshold for read quality filtering, default=10
  -d MINDIST, --mindist MINDIST
                        Minimum distance (in base pairs), default=5000
  -r RESOLUTIONS [RESOLUTIONS ...], --resolutions RESOLUTIONS [RESOLUTIONS ...]
                        Default resolutions for the matrices default=300 150 50
  --keep-intermediate-files
                        Keep intermediate files from each stage

Minimal usage example: te_hic -a read1.bam -b read2.bam -g genome
```

## HiC example

A toy example is included in te_hic/test/HiC-analysis (the 'HiC-analysis % ' is the user prompt, type
in the commands after the % to get it it run)

```
HiC-analysis % te_hic -1 SRR1030718.2_1.100k.bam -2 SRR1030718.2_2.100k.bam -g hg38 
INFO    : Arguments:
INFO    :   Read1: SRR1030718.2_1.100k.bam
INFO    :   Read2: SRR1030718.2_2.100k.bam
INFO    :   Genome: hg38
INFO    :   Label: te_hic (default is te_hic)
INFO    :   Quality thresold: 10 (default is 10)
INFO    :   Minimum contact distance: 5000 (default is 5000)
INFO    :   Matrix resolutions to build: [300, 150, 50] kbp (default is [300, 150, 50])
INFO    : valid chromosome are: chr1, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr2, chr20, chr21, chr22, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chrX, chrY
INFO    : Loaded '/Users/andrew/Tools/te_hic/tehiclib/../genome/hg38_te_genome_freqs.glb' binary file with 1058 items
INFO    : Stage 1: Collect valid read pairs
INFO    : Started SRR1030718.2_1.100k.bam and SRR1030718.2_2.100k.bam
[E::idx_find_and_load] Could not retrieve index file for 'SRR1030718.2_1.100k.bam'
[E::idx_find_and_load] Could not retrieve index file for 'SRR1030718.2_2.100k.bam'
INFO    : Processed: 99,970 reads in total, removed 47,357 (47.4%) reads by preduplicate removal
INFO    : Sorting temp file
INFO    : Stage 1 stats:
INFO    :   Aligned:
INFO    :     Reads processed           : 99,970
INFO    :     Correctly paired          : 85,809 (85.83%)
INFO    :   Rejected reads:
INFO    :     Low quality               : 4,213 (4.21%)
INFO    :     One pair aligned          : 12,252 (12.26%)
INFO    :     Not canonical chromosome  : 312 (0.31%)
INFO    :     No pairs aligned          : 1,909 (1.91%)
INFO    :   Rejected reads (by criteria):
INFO    :     Too close                 : 28,646 (28.65%)
INFO    :   Final:
INFO    :     Kept reads                : 52,613 (52.63%)
INFO    :     Trans links               : 19,133 (19.14%)
INFO    :     [Note that these numbers below include PCR duplicates]
INFO    :     Kept short-range (<20kb)  : 2,418 (2.42%)
INFO    :     Kept long-range (>20kb)   : 31,087 (31.10%)
INFO    :     [Sorted, unique after full duplicate removal]
INFO    :     Kept long-range (>20kb)   : 52,613 (52.63%)
INFO    : Took 0.8s
INFO    : Stage 2: Assign to hg38 genome features
INFO    : Loaded '/Users/andrew/Tools/te_hic/tehiclib/../genome/hg38_tes.glb' binary file with 5665930 items
INFO    : Processed: 0
INFO    : Processed: 52,612 reads in total
INFO    : Saving BEDPE files
INFO    : Took 29.4s
INFO    : Stage 3: Quantify links
INFO    : Measures anchors...
INFO    :   TE <-> TE : 15,578 (29.61%)
INFO    :   TE <-> -- : 25,300 (48.09%)
INFO    :   -- <-> -- : 11,735 (22.30%)
INFO    : Saved 'stage3.te_hic_te-nn_anchor_frequencies.tsv'
INFO    : Took 1.0s
INFO    : Stage 4: Build Matrices
INFO    : Building in-memory matrices for resolution 300000 bp
INFO    : Saved BED bins: "matrices_te_hic/300000/te_hic_300000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/300000/te_hic_300000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/300000/te_hic_300000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/300000/te_hic_300000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/300000/te_hic_300000.nnnn.raw.matrix"
INFO    : Building in-memory matrices for resolution 150000 bp
INFO    : Saved BED bins: "matrices_te_hic/150000/te_hic_150000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/150000/te_hic_150000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/150000/te_hic_150000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/150000/te_hic_150000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/150000/te_hic_150000.nnnn.raw.matrix"
INFO    : Building in-memory matrices for resolution 50000 bp
INFO    : Saved BED bins: "matrices_te_hic/50000/te_hic_50000_abs.bed"
INFO    : Saved All matrix: "matrices_te_hic/50000/te_hic_50000.all.raw.matrix"
INFO    : Saved TE <=> TE matrix: "matrices_te_hic/50000/te_hic_50000.tete.raw.matrix"
INFO    : Saved TE <=> non-TE matrix: "matrices_te_hic/50000/te_hic_50000.tenn.raw.matrix"
INFO    : Saved non-TE <=> non-TE matrix: "matrices_te_hic/50000/te_hic_50000.nnnn.raw.matrix"
INFO    : Took 1.0s
INFO    : In total took 32.2s

HiC-analysis % ls -l
total 16960
drwxr-xr-x  5 andrew  staff      160 Dec 21  2022 matrices_te_hic
-rw-r--r--  1 andrew  staff   281469 Oct 17 12:27 stage1.te_hic.5000.frags.bed.gz
-rw-r--r--@ 1 andrew  staff  1234822 Oct 17 12:27 stage2.all.te_hic.bedpe.gz
-rw-r--r--  1 andrew  staff   232363 Oct 17 12:27 stage2.nnnn.te_hic.bedpe.gz
-rw-r--r--  1 andrew  staff   589225 Oct 17 12:27 stage2.tenn.te_hic.bedpe.gz
-rw-r--r--  1 andrew  staff   401417 Oct 17 12:27 stage2.tete.te_hic.bedpe.gz
-rw-r--r--  1 andrew  staff       93 Oct 17 12:27 stage3.te_hic_crude_measures.txt
-rw-r--r--  1 andrew  staff    79178 Oct 17 12:27 stage3.te_hic_te-nn_anchor_frequencies.tsv
-rw-r--r--  1 andrew  staff  1357255 Oct 17 12:27 stage3.te_hic_te-te_anchor_frequencies.tsv


HiC-analysis % ll matrices_te_hic 
total 0
drwxr-xr-x  7 andrew  staff  224 Dec 21 08:31 150000
drwxr-xr-x  7 andrew  staff  224 Dec 21 08:31 300000
drwxr-xr-x  7 andrew  staff  224 Dec 21 08:31 50000

HiC-analysis % ll matrices_te_hic/150000 
total 3520
-rw-r--r--  1 andrew  staff  585827 Dec 21 08:31 te_hic_150000.all.raw.matrix
-rw-r--r--  1 andrew  staff  130043 Dec 21 08:31 te_hic_150000.nnnn.raw.matrix
-rw-r--r--  1 andrew  staff  287113 Dec 21 08:31 te_hic_150000.tenn.raw.matrix
-rw-r--r--  1 andrew  staff  181391 Dec 21 08:31 te_hic_150000.tete.raw.matrix
-rw-r--r--  1 andrew  staff  603430 Dec 21 08:31 te_hic_150000_abs.bed

HiC-analysis % wc *.tsv
     904    5428   79178 stage3.te_hic_te-nn_anchor_frequencies.tsv
   11570   69423 1357255 stage3.te_hic_te-te_anchor_frequencies.tsv
   12474   74851 1436433 total
   
HiC-analysis % head stage3.te_hic_te-nn_anchor_frequencies.tsv 
name	count	genome_count	genome_percent	RPM	RPM per kbp of TE
DNA:DNA:Eulor1	1	11998	0.0003874509893470593	19.007793195210038	1.5842468074020701
DNA:DNA:Eulor9C	1	9179	0.0002964171221217417	19.007793195210038	2.070791283931805
DNA:DNA:MER125	2	15748	0.0005085496066208942	38.015586390420076	2.4139945637808022
DNA:DNA:MER126	3	41040	0.00132530326744485	57.0233795856301	1.3894585669013184
DNA:DNA:MER135	5	78626	0.0025390666351393465	95.03896597605019	1.2087473097455064
DNA:DNA:MER136	1	4736	0.000152939480375702	19.007793195210038	4.013469846961579
DNA:MULE-MuDR:Ricksha	3	119602	0.00386230315284939	57.0233795856301	0.47677613740263625
DNA:MULE-MuDR:Ricksha_0	1	190682	0.006157687077069174	19.007793195210038	0.09968320657015364
DNA:MULE-MuDR:Ricksha_a	1	24538	0.0007924047655107635	19.007793195210038	0.7746268316574308

HiC-analysis % head stage3.te_hic_te-te_anchor_frequencies.tsv
TE1	TE2	RPM	RPM per kbp TE	TE1_genome_freq	TE2_genome_freq
DNA:DNA:Eulor1	LINE:L1:L1MD1	19.007793195210038	0.004780921880629596	3.874509893470593e-06	0.0012800159368105426
DNA:DNA:MER126	DNA:TcMar-Tigger:Tigger9a	19.007793195210038	0.08175957569213382	1.32530326744485e-05	6.182294316099218e-05
DNA:DNA:MER126	LINE:CR1:L3	19.007793195210038	0.0021427946949933996	1.32530326744485e-05	0.002851314414370416
DNA:DNA:MER126	LINE:L2:L2	19.007793195210038	0.0012459431318053477	1.32530326744485e-05	0.0049132799464707685
DNA:DNA:MER126	LTR:ERVL-MaLR:MLT1K	19.007793195210038	0.004841691471064424	1.32530326744485e-05	0.001254522901761347
DNA:DNA:MER130	LTR:ERV1:LTR29	19.007793195210038	0.07713074873480352	1.0474223070071568e-05	6.91072671872479e-05
DNA:DNA:MER135	LINE:L1:L1ME3G	19.007793195210038	0.004568449665320416	2.5390666351393466e-05	0.0013182117324172945
DNA:DNA:MER135	SINE:Alu:AluJo	19.007793195210038	0.0010892240652607514	2.5390666351393466e-05	0.005609978698636967
DNA:DNA:MER135	SINE:Alu:AluSc5	19.007793195210038	0.00985529706993046	2.5390666351393466e-05	0.0005974398668556419
```

## License

Release License:
MIT license:
Copyright (C) 2019-2024 Andrew Hutchins
    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
Except as contained in this notice, the name(s) of the above copyright holders shall not be used in advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization.
    
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


