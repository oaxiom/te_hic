
getloops -m bed_to_gene -g hg38 -i testbedpe.chr10.bedpe.gz -p hg38_random.chr10.bed.gz -o testout.tsv.gz

# Should give:

INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['hg38_random.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "['hg38']"
INFO    :   mode: "bed_to_gene"
INFO    :   window: 5000
INFO    : Loaded '/Users/andrew/Tools/te_hic/bin/../genome/hg38_glb_gencode_promoters.glb' binary file with 1160903 items
INFO    : Loaded /Users/andrew/Tools/te_hic/bin/../genome/hg38_glb_gencode_promoters.glb
INFO    : Found 432 BED peaks
INFO    : Found 1,160,902 transcripts
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   5838 (92.6%) loops have 1 reads
INFO    :   393 (6.2%) loops have 2 reads
INFO    :   51 (0.8%) loops have 3 reads
INFO    :   13 (0.2%) loops have 4 reads
INFO    :   4 (0.1%) loops have 5 reads
INFO    :   1 (0.0%) loops have 6 reads
INFO    :   4 (0.1%) loops have 7 reads
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

getloops -m bed_to_bed -i testbedpe.chr10.bedpe.gz -p CTCF_constitutive.chr10.bed.gz -o testout.tsv.gz

# Should give:
getloops -m bed_to_gene -g hg38 -i testbedpe.chr10.bedpe.gz -p CTCF_constitutive.chr10.bed.gz -o testout.tsv.gz 
/Users/andrew/Library/Python/3.8/lib/python/site-packages/umap/spectral.py:4: NumbaDeprecationWarning: No direct replacement for 'numba.targets' available. Visit https://gitter.im/numba/numba-dev to request help. Thanks!
  import numba.targets
INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['CTCF_constitutive.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "['hg38']"
INFO    :   mode: "bed_to_gene"
INFO    :   window: 5000
INFO    : Loaded '/Users/andrew/Tools/te_hic/bin/../genome/hg38_glb_gencode_promoters.glb' binary file with 1160903 items
INFO    : Loaded /Users/andrew/Tools/te_hic/bin/../genome/hg38_glb_gencode_promoters.glb
INFO    : Found 1,263 BED peaks
INFO    : Found 1,160,902 transcripts
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   22836 (90.2%) loops have 1 reads
INFO    :   1907 (7.5%) loops have 2 reads
INFO    :   363 (1.4%) loops have 3 reads
INFO    :   101 (0.4%) loops have 4 reads
INFO    :   62 (0.2%) loops have 5 reads
INFO    :   21 (0.1%) loops have 6 reads
INFO    :   9 (0.0%) loops have 7 reads
INFO    :   6 (0.0%) loops have 8 reads
INFO    :   2 (0.0%) loops have 9 reads
INFO    :   3 (0.0%) loops have 10 reads
INFO    :   0 (0.0%) loops have 11 reads
INFO    :   0 (0.0%) loops have 12 reads
INFO    :   2 (0.0%) loops have 13 reads
INFO    :   0 (0.0%) loops have 14 reads
INFO    :   0 (0.0%) loops have 15 reads
INFO    :   0 (0.0%) loops have 16 reads
INFO    :   0 (0.0%) loops have 17 reads
INFO    :   0 (0.0%) loops have 18 reads
INFO    :   0 (0.0%) loops have 19 reads
INFO    :   0 (0.0%) loops have 20+ reads
INFO    : Saved testout.tsv.gz
