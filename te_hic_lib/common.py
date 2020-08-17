
import sys

valid_assemblies = {'mm10', 'hg38'}
valid_modes = {'bed_to_bed', 'bed_to_genes', 'genes_to_tes', 'bed_to_tes', 'bed_to_tes_to_genes'}

def print_genomes(log=None):
    log.info('  Valid Genome assemblies are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

