
import sys

valid_assemblies = {'hg38'} # mm10 is temporarily disabled;
valid_modes = {'bed_to_bed', 'bed_to_gene', 'genes_to_te', 'bed_to_te', 'bed_to_tesandgenes'}

def print_genomes(log=None):
    log.info('  Valid Genome assemblies are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

