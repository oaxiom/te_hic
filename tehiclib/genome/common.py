'''

Package for genome annotations

'''

import os, gzip

valid_assemblies = {
    'hg38',
    'mm10',
    'mm39',
    'rn7',
    'xenTro10',
    }

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
    # hg38, mm10

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

genome_options = {
    'hg38': {
        'download': "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
        'chrom_cleaner': clean_chroms_animal,
        },
    'mm10': {
        'download': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz',
        'chrom_cleaner': clean_chroms_animal,
        },
    'mm39': {
        'download': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz',
        'chrom_cleaner': clean_chroms_animal,
        },
    'rn7': {
        'download': 'https://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz',
        'chrom_cleaner': clean_chroms_animal,
        }
    }


