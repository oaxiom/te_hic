#!/usr/bin/env python3

'''

Build the combined gencode and repeat data for te_hic

This is for the mm10 genome

'''

import sys
sys.path.append('../')
from tehiclib import delayedlist, progressbar, genelist

gtf = {
    #"feature_type": 1,
    "feature": 2,
    "gtf_decorators": 8,
    "commentlines": "#",
    "loc": "location(chr=column[0], left=column[3], right=column[4])",
    "strand": 6,
    "skiplines": -1,
    "force_tsv": True
    }

rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
    'repName': 10, 'repClass': 11, 'repFamily': 12}

chr_set = frozenset(['X', 'Y'] + [str(i) for i in range(1, 22)])

repeats = delayedlist(filename='hg38_rmsk.txt.gz', gzip=True, format=rmsk_track_form)
gencode = delayedlist('gencode.v29.annotation.gtf.gz', gzip=True, format=gtf)

keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon'])

added = 0

newl = []
promoters = []

print('Repeats')
p = progressbar(len(repeats))
for idx, item in enumerate(repeats):
    if item['repClass'] not in keep_classes:
        continue
    if str(item['loc']['chr']) not in chr_set:
        continue

    name = f"{item['repClass']}:{item['repFamily']}:{item['repName']}"

    newentry = {'loc': item['loc'],
        'name': name,
        'type': 'TE',
        'ensg': name
        }
    newl.append(newentry)

    added += 1

    #if idx > 100000:
    #    break
    p.update(idx)

print(f'\nAdded {added:,} features')

keep_gene_types = set(('protein_coding', 'lincRNA', 'lncRNA'))

print('Gencode')
p = progressbar(len(gencode))
for idx, item in enumerate(gencode):
    if item['feature'] != 'exon':
        continue
    if item['gene_type'] not in keep_gene_types:
        continue
    if item['loc'].loc['chr'] not in chr_set:
        continue

    newentry = {
        'loc': item['loc'],
        'name': item['gene_name'],
        'type': item['gene_type'],
        'ensg': item['gene_id'].split('.')[0],
        }
    newl.append(newentry)
    added += 1

    if item['strand'] == '+':
        prom_locs = {
            'loc': item['loc'].pointLeft(),
            'name': item['gene_name'],
            'type': item['gene_type'],
            'ensg': newentry['ensg'],
            'enst': item['transcript_id'].split('.')[0],
            }
    elif item['strand'] == '-':
        prom_locs = {
            'loc': item['loc'].pointRight(),
            'name': item['gene_name'],
            'type': item['gene_type'],
            'ensg': newentry['ensg'],
            'enst': item['transcript_id'].split('.')[0],
            }
    else:
        1/0

    promoters.append(prom_locs)

    #if idx > 100000:
    #    break

    p.update(idx)

print(f'\nAdded {added:,} features')

gl = genelist()
gl.load_list(promoters)
gl.save('hg38_glb_gencode_promoters.glb')

gl = genelist()
gl.load_list(newl)
gl.save('hg38_glb_gencode_tes.glb')
genome = gl

hg38_genome_size = 3096649726 # http://asia.ensembl.org/Homo_sapiens/Info/Annotation

tes = {}

for idx, te in enumerate(genome):
    if ':' not in te['name']:
        continue # omit the genes
    if '?' in te:
        continue

    if te['name'] not in tes:
        tes[te['name']] = 0
    tes[te['name']] += len(te['loc'])

newl = []
for k in tes:
    newe = {'name': k,
        'genome_count': tes[k],
        'genome_percent': tes[k] / hg38_genome_size * 100.0}
    newl.append(newe)

gl = genelist()
gl.load_list(newl)
gl.sort('name')
gl.saveTSV('hg38_te_genome_freqs.tsv', key_order=['name', 'genome_count', 'genome_percent'])
gl.save('hg38_te_genome_freqs.glb')



