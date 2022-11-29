
import sys

valid_assemblies = {'hg38'} # mm10 is temporarily disabled;
valid_modes = {'bed_to_bed', 'bed_to_gene', 'genes_to_te', 'bed_to_te', 'bed_to_tesandgenes'}

def print_genomes(log=None):
    log.info('  Valid Genome assemblies are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

# http://asia.ensembl.org/Homo_sapiens/Info/Annotation
hg38_genome_size      = 3096649726
mm10_genome_size      = 2728222451
mm39_genome_size      = 2728222451
danRer11_genome_size  = 1373471384 # GRCz11
dm6_genome_size       =  143726002 # BDGP6_32
rn7_genome_size       = 2647915728 # RGSC 6.0
xenTro10_genome_size  = 1451301209 # UCB_Xtro_10
