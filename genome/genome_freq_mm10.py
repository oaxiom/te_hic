#!/usr/bin/env python3
'''

measure the % and number ofbase pairs each type of TE occupuies in the genome

'''

from glbase3 import *

genome = glload('mm10_glb_gencode_tes.glb')

mm10_genome_size = 2730871774 # http://asia.ensembl.org/Mus_musculus/Info/Annotation

tes = {}

for te in genome:
    if '?' in te:
        continue
    if ':' not in te['name']:
        continue # omit the genes
    if te['name'] not in tes:
        tes[te['name']] = 0
    tes[te['name']] += len(te['loc'])

newl = []
for k in tes:
    newe = {'name': k,
        'count': tes[k],
        'perc': tes[k] / mm10_genome_size * 100.0}
    newl.append(newe)

gl = genelist()
gl.load_list(newl)
gl.sort('name')
gl.saveTSV('mm10_te_genome_freqs.tsv', key_order=['name', 'count', 'perc'])
gl.save('mm10_te_genome_freqs.glb')


