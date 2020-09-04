
getloops -m bed_to_gene -g mm10 -i testbedpe.chr10.bedpe.gz -p hg38_random.chr10.bed.gz -o testout.tsv.gz

# Should give:

INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['hg38_random.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "['mm10']"
INFO    :   mode: "bed_to_gene"
INFO    :   window: 5000
INFO    : Loaded '/Users/andrew/Tools/te_hic/bin/../genome/mm10_glb_gencode_promoters.glb' binary file with 763710 items
INFO    : Loaded /Users/andrew/Tools/te_hic/bin/../genome/mm10_glb_gencode_promoters.glb
INFO    : Found 432 BED peaks
INFO    : Shrunk peaks from 432 to 334 by removing duplicate bins and common bins
INFO    : Shrunk genes from 763,709 to 117,675 by removing duplicate bins and common bins
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   2944 (92.5%) loops have 1 reads
INFO    :   194 (6.1%) loops have 2 reads
INFO    :   35 (1.1%) loops have 3 reads
INFO    :   6 (0.2%) loops have 4 reads
INFO    :   2 (0.1%) loops have 5 reads
INFO    :   2 (0.1%) loops have 6 reads
INFO    :   0 (0.0%) loops have 7 reads
INFO    :   1 (0.0%) loops have 8 reads
INFO    :   0 (0.0%) loops have 9 reads
INFO    :   0 (0.0%) loops have 10 reads
INFO    :   0 (0.0%) loops have 11 reads
INFO    :   0 (0.0%) loops have 12 reads
INFO    :   0 (0.0%) loops have 13 reads
INFO    :   0 (0.0%) loops have 14 reads
INFO    :   0 (0.0%) loops have 15 reads
INFO    :   0 (0.0%) loops have 16 reads
INFO    :   0 (0.0%) loops have 17 reads
INFO    :   0 (0.0%) loops have 18 reads
INFO    :   0 (0.0%) loops have 19 reads
INFO    :   0 (0.0%) loops have 20+ reads
INFO    : Saved testout.tsv.gz

getloops -m bed_to_gene -g mm10 -i testbedpe.chr10.bedpe.gz -p CTCF_constitutive.chr10.bed.gz -o testout.tsv.gz 

INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['CTCF_constitutive.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "['mm10']"
INFO    :   mode: "bed_to_gene"
INFO    :   window: 5000
INFO    : Loaded '/Users/andrew/Tools/te_hic/bin/../genome/mm10_glb_gencode_promoters.glb' binary file with 763710 items
INFO    : Loaded /Users/andrew/Tools/te_hic/bin/../genome/mm10_glb_gencode_promoters.glb
INFO    : Found 1,263 BED peaks
INFO    : Shrunk peaks from 1,263 to 720 by removing duplicate bins and common bins
INFO    : Shrunk genes from 763,709 to 117,242 by removing duplicate bins and common bins
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   10367 (89.5%) loops have 1 reads
INFO    :   929 (8.0%) loops have 2 reads
INFO    :   183 (1.6%) loops have 3 reads
INFO    :   53 (0.5%) loops have 4 reads
INFO    :   24 (0.2%) loops have 5 reads
INFO    :   15 (0.1%) loops have 6 reads
INFO    :   3 (0.0%) loops have 7 reads
INFO    :   1 (0.0%) loops have 8 reads
INFO    :   0 (0.0%) loops have 9 reads
INFO    :   1 (0.0%) loops have 10 reads
INFO    :   1 (0.0%) loops have 11 reads
INFO    :   1 (0.0%) loops have 12 reads
INFO    :   1 (0.0%) loops have 13 reads
INFO    :   0 (0.0%) loops have 14 reads
INFO    :   0 (0.0%) loops have 15 reads
INFO    :   0 (0.0%) loops have 16 reads
INFO    :   0 (0.0%) loops have 17 reads
INFO    :   0 (0.0%) loops have 18 reads
INFO    :   0 (0.0%) loops have 19 reads
INFO    :   0 (0.0%) loops have 20+ reads
INFO    : Saved testout.tsv.gz