'''

Package for genome annotations

'''

import os

valid_assemblies = {'mm10', 'hg38'}

from .make import make_index

def check_genome_done(genome):
    script_path = os.path.dirname(os.path.realpath(__file__))
    # Check all three genome idxs are avaialable;

    if os.path.exists(f'{script_path}/../../genome/{genome}_tes.glb'):
        #if os.path.exists(f'{script_path}/../../genome/{genome}_promoters.glb'):
        if os.path.exists(f'{script_path}/../../genome/{genome}_te_genome_freqs.glb'):
            return True

    return False