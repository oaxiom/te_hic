'''

Package for genome annotations

'''

import os

from .make import make_index
from .common import valid_assemblies

def check_genome_done(genome):
    script_path = os.path.dirname(os.path.realpath(__file__))
    # Check all genome idxs are avaialable;

    if os.path.exists(f'{script_path}/../../genome/{genome}_tes.glb'):
        if os.path.exists(f'{script_path}/../../genome/{genome}_te_genome_freqs.glb'):
            return True

    return False
