
getContacts -m bed_to_bed -i testbedpe.chr10.bedpe.gz -p hg38_random.chr10.bed.gz -o testout.tsv.gz

# Should give:
INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['hg38_random.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "None"
INFO    :   mode: "bed_to_bed"
INFO    :   window: 2000
INFO    : Found 432 BED peaks
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   9 (100.0%) loops have 1 reads
INFO    :   0 (0.0%) loops have 2 reads
INFO    :   0 (0.0%) loops have 3 reads
INFO    :   0 (0.0%) loops have 4 reads
INFO    :   0 (0.0%) loops have 5 reads
INFO    :   0 (0.0%) loops have 6 reads
INFO    :   0 (0.0%) loops have 7 reads
INFO    :   0 (0.0%) loops have 8 reads
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

getContacts -m bed_to_bed -i testbedpe.chr10.bedpe.gz -p CTCF_constitutive.chr10.bed.gz -o testout.tsv.gz


# Should give:
INFO    : getloops
INFO    : Arguments:
INFO    :   inreadsbedpe: ['testbedpe.chr10.bedpe.gz']
INFO    :   inpeaksbed: ['CTCF_constitutive.chr10.bed.gz']
INFO    :   outtsv: ['testout.tsv.gz']
INFO    :   genome: "None"
INFO    :   mode: "bed_to_bed"
INFO    :   window: 2000
INFO    : Found 1,263 BED peaks
INFO    : 1,000,000 reads processed
INFO    : Histogram of loops:
INFO    :   666 (89.3%) loops have 1 reads
INFO    :   72 (9.7%) loops have 2 reads
INFO    :   6 (0.8%) loops have 3 reads
INFO    :   2 (0.3%) loops have 4 reads
INFO    :   0 (0.0%) loops have 5 reads
INFO    :   0 (0.0%) loops have 6 reads
INFO    :   0 (0.0%) loops have 7 reads
INFO    :   0 (0.0%) loops have 8 reads
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
