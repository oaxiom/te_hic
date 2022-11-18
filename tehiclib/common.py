
import sys

valid_assemblies = {'mm10', 'hg38'}
valid_modes = {'bed_to_bed', 'bed_to_gene', 'genes_to_te', 'bed_to_te', 'bed_to_tesandgenes'}

def print_genomes(log=None):
    log.info('  Valid Genome assemblies are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

