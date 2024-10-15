#!/usr/bin/env python3

'''

Build the combined gencode and repeat data for te_hic

This is for the mm10 genome

'''

import sys, subprocess, os, gzip
sys.path.append('../')
from ..miniglbase3 import delayedlist, progressbar, genelist
from .common import genome_sizes, clean_chroms_animal

# IS it possible to avoid hardcoding these?
# Also can include options for specific genome filtering?
genome_options = {
    'hg38': {
        'download': "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
        'chrom_cleaner': clean_chroms_animal,
        },
    'mm10': {
        'download': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz',
        'chrom_cleaner': clean_chroms_animal,
        }
    }

def make_index(genome, log):
    script_path = os.path.dirname(os.path.realpath(__file__))

    # Download data:
    chrom_sizes_path = f'{script_path}/../../genome/{genome}.chromSizes.gz'
    rmsk_path = f'{script_path}/../../genome/{genome}.rmsk.txt.gz'
    #annotation_path = f'{script_path}/../../genome/{genome}.annotation.txt.gz'

    subprocess.run(f'wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/chromInfo.txt.gz -O {chrom_sizes_path}', shell=True)
    subprocess.run(f'wget -c http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz -O {rmsk_path}', shell=True)
    #subprocess.run(f"wget -c {genome_options[genome]['download']} -O {annotation_path}", shell=True)

    if genome_options[genome]['chrom_cleaner']:
        genome_options[genome]['chrom_cleaner'](genome)

    rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
        'repName': 10, 'repClass': 11, 'repFamily': 12}

    # TODO: Fix for other species?
    chr_set = frozenset(['X', 'Y'] + [str(i) for i in range(1, 23)])

    ###### Repeats table;
    repeats = delayedlist(filename=rmsk_path, gzip=True, format=rmsk_track_form)

    # TODO: Needs to be expanded for other species?
    keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon'])

    added = 0

    newl = []
    promoters = []

    log.info('Adding repeat annotations')
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

        p.update(idx)

    log.info(f'\nAdded {added:,} features')
    del repeats

    ###### Annotation table
    gtf = {
        #"feature_type": 1,
        "feature": 2,
        "gtf_decorators": 8,
        "commentlines": "#",
        "loc": "location(chr=column[0], left=column[3], right=column[4])",
        "strand": 6,
        "skiplines": -1,
        "force_tsv": True,
        }

    gencode = delayedlist(annotation_path, gzip=True, format=gtf)
    keep_gene_types = set(('protein_coding', 'lincRNA', 'lncRNA'))


    log.info('Adding gene annotations')
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
            1/0 # Panic!

        promoters.append(prom_locs)

        p.update(idx)

    log.info(f'\nAdded {added:,} features')

    gl = genelist()
    gl.load_list(promoters)
    gl.save(f'{script_path}/../../genome/{genome}_promoters.glb')

    gl_tes = genelist()
    gl_tes.load_list(newl)
    gl_tes.save(f'{script_path}/../../genome/{genome}_tes.glb') # Includes genes and TEs;

    log.info(f'Genome annotation has {len(gl_tes):,} features in total')

    # Count all the TE types
    tes = {}

    for idx, te in enumerate(gl_tes):
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
            'genome_percent': tes[k] / genome_sizes[genome] * 100.0}
        newl.append(newe)

    gl = genelist()
    gl.load_list(newl)
    gl.sort('name')
    gl.saveTSV(f'{script_path}/../../genome/{genome}_te_genome_freqs.tsv', key_order=['name', 'genome_count', 'genome_percent'])
    gl.save(f'{script_path}/../../genome/{genome}_te_genome_freqs.glb')



