'''

Package for genome annotations

'''

import os

valid_assemblies = {'mm10', 'hg38'}

from .make import make_index

def check_genome_done(genome):
    # Check all three genome idxs are avaialable;

    if not os.path.exists(f'../../genome/{genome}_glb_gencode_tes.glb'):
        return False

    if not os.path.exists(f'../../genome/{genome}_glb_gencode_promoters.glb'):
        return False

    if not os.path.exists(f'../../genome/{genome}_te_genome_freqs.glb'):
        return False
