#!/usr/bin/env python3

'''

Build the combined gencode and repeat data for te_hic

This is for the mm10 genome

'''

import sys
sys.path.append('../tehiclib')
from miniglbase3 import delayedlist, progressbar

gtf = {
    #"feature_type": 1,
    "feature": 2,
    "gtf_decorators": 8,
    "commentlines": "#",
    "loc": "location(chr=column[0], left=column[3], right=column[4])",
    "strand": 6,
    "skiplines": -1,
    "force_tsv": True}

rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
    'repName': 10, 'repClass': 11, 'repFamily': 12}

chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

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

    newentry = {'loc': item['loc'],
        'name': '%s:%s:%s' % (item['repName'], item['repFamily'], item['repClass']),
        'type': 'TE',
        'ensg': '%s:%s:%s' % (item['repName'], item['repFamily'], item['repClass'])
        }
    newl.append(newentry)
    #print(newentry)
    added += 1

    #if idx > 100000:
    #    break
    p.update(idx)
print('\nAdded %s features' % added)

print('Gencode')
p = progressbar(len(gencode))
for idx, item in enumerate(gencode):
    if item['feature'] != 'exon':
        continue
    if item['gene_type'] not in ('protein_coding', 'lincRNA'):
        continue
    if item['loc']['chr'] not in chr_set:
        continue

    newentry = {'loc': item['loc'],
        'name': item['gene_name'],
        'type': item['gene_type'],
        'ensg': item['gene_id'].split('.')[0],
        }
    newl.append(newentry)
    added += 1

    if item['strand'] == '+':
        prom_locs = {'loc': item['loc'].pointLeft(),
            'name': item['gene_name'],
            'type': item['gene_type'],
            'ensg': item['gene_id'].split('.')[0],
            'enst': item['transcript_id'].split('.')[0],
            }
    elif item['strand'] == '-':
        prom_locs = {'loc': item['loc'].pointRight(),
            'name': item['gene_name'],
            'type': item['gene_type'],
            'ensg': item['gene_id'].split('.')[0],
            'enst': item['transcript_id'].split('.')[0],
            }
    else:
        1/0

    promoters.append(prom_locs)

    #if idx > 100000:
    #    break

    p.update(idx)

print('\nAdded %s features' % added)

gl = genelist()
gl.load_list(newl)
gl.save('hg38_glb_gencode_tes.glb')

gl = genelist()
gl.load_list(promoters)
gl.save('hg38_glb_gencode_promoters.glb')

print(gl)
