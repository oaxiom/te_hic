'''

Package for genome annotations

'''

import os, gzip

genome_sizes = dict(
    hg38      = 3096649726,
    mm10      = 2728222451,
    mm39      = 2728222451,
    danRer11  = 1373471384, # GRCz11
    dm6       =  143726002, # BDGP6_32
    rn7       = 2647915728, # RGSC 6.0
    xenTro10  = 1451301209, # UCB_Xtro_10
    )

# Cleaners for chromsome names to define canonical;
def clean_chroms_animal(genome):
    # Confirmed to work for:
    # hg38
    
    script_path = os.path.dirname(os.path.realpath(__file__))
    
    out = open(f'{script_path}/../../genome/{genome}.chromSizes.clean', 'wt')
    
    with gzip.open(f'{script_path}/../../genome/{genome}.chromSizes.gz', 'rt') as inp:
        for chrom in inp:
            if '_alt' in chrom: continue
            if '_fix' in chrom: continue
            if '_random' in chrom: continue
            if '_KI' in chrom: continue
            if chrom.startswith('chrUn'): continue
            if chrom.startswith('chrM'): continue
            
            out.write(f'{chrom}')
    out.close()