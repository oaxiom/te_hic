#!/usr/bin/env python3

'''

Build the combined gencode and repeat data for te_hic

This is for the mm10 genome

'''

import sys, subprocess, os, gzip
sys.path.append('../')
from ..miniglbase3 import delayedlist, progressbar, genelist
from .common import genome_sizes, clean_chroms_animal, genome_options

# IS it possible to avoid hardcoding these?
# Also can include options for specific genome filtering?


def make_index(genome, log):
    script_path = os.path.dirname(os.path.realpath(__file__))

    # Download data:
    chrom_sizes_path = f'{script_path}/../../genome/{genome}.chromSizes.gz'
    rmsk_path = f'{script_path}/../../genome/{genome}.rmsk.txt.gz'
    annotation_path = f'{script_path}/../../genome/{genome}.annotation.txt.gz'

    subprocess.run(f'wget -c http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/chromInfo.txt.gz -O {chrom_sizes_path}', shell=True)
    subprocess.run(f'wget -c http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz -O {rmsk_path}', shell=True)
    subprocess.run(f"wget -c {genome_options[genome]['download']} -O {annotation_path}", shell=True)

    if genome_options[genome]['chrom_cleaner']:
        valid_chroms = genome_options[genome]['chrom_cleaner'](genome)

    rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
        'repName': 10, 'repClass': 11, 'repFamily': 12}

    ###### Repeats table;
    repeats = delayedlist(filename=rmsk_path, gzip=True, format=rmsk_track_form)

    # TODO: Needs to be expanded for other species?
    keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'RNA', 'Retroposon'])
    classes_seen = {}
    chromosomes_seen_but_not_used = set()

    added = 0

    newl = []
    promoters = []

    log.info('Adding repeat annotations')
    p = progressbar(len(repeats))
    for idx, item in enumerate(repeats):
        if item['repClass'] not in keep_classes:
            if '?' not in item['repClass']:
                if item['repClass'] not in classes_seen: 
                    classes_seen[item['repClass']] = 0
                classes_seen[item['repClass']] += 1
            continue

        if str(item['loc'].chrom) not in valid_chroms:
            chromosomes_seen_but_not_used.add(item['loc'].chrom)
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

    print() # tidy up after progress bar
    log.info('Chromsomes seen but not used (Scaffolds with "_" are not included):')
    log.info(', '.join([c for c in chromosomes_seen_but_not_used if '_' not in c]))
    
    log.info('TE types seen but not added:')
    for seen_te in classes_seen:
        log.info(f'    {seen_te} ({classes_seen[seen_te]:,} TE loci)')

    log.info(f'Added {added:,} features')
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
    keep_gene_types = set(('protein_coding', 'lincRNA', 'lncRNA', 'ncRNA'))
    gene_types_seen = {}

    key_gene_type = None
    key_gene_name = None

    log.info('Adding gene annotations')
    p = progressbar(len(gencode))
    for idx, item in enumerate(gencode):
        if item['feature'] != 'exon':
            continue
        if item['loc'].chrom not in valid_chroms:
            continue

        # guess the gene_biotype and name;
        if not key_gene_type:
            if 'gene_biotype' in item: key_gene_type = 'gene_biotype'
            if 'gene_type' in item: key_gene_type = 'gene_type'
        if not key_gene_name:
            if 'gene_name' in item: key_gene_name = 'gene_name'

        if item[key_gene_type] not in keep_gene_types:
            if item[key_gene_type] not in gene_types_seen:
                gene_types_seen[item[key_gene_type]] = 1  
            gene_types_seen[item[key_gene_type]] += 1
            continue

        # type is always valid (right?) but gene name might not be.
        gene_type = item[key_gene_type]
        try:
            gene_name = item[key_gene_name]
        except KeyError:
            gene_name = item['gene_id'].split('.')[0]

        newentry = {
            'loc': item['loc'],
            'name': gene_name,
            'type': gene_type,
            'ensg': item['gene_id'].split('.')[0],
            }
        newl.append(newentry)
        added += 1

        if item['strand'] == '+':
            prom_locs = {
                'loc': item['loc'].pointLeft(),
                'name': gene_name,
                'type': gene_type,
                'ensg': newentry['ensg'],
                'enst': item['transcript_id'].split('.')[0],
                }
        elif item['strand'] == '-':
            prom_locs = {
                'loc': item['loc'].pointRight(),
                'name': gene_name,
                'type': gene_type,
                'ensg': newentry['ensg'],
                'enst': item['transcript_id'].split('.')[0],
                }
        else:
            1/0 # Panic!

        promoters.append(prom_locs)

        p.update(idx)

    print() # tidy up after progressbar
    log.info('gene_types seen but not added:')
    for seen_gene in gene_types_seen:
        log.info(f'    {seen_gene} ({gene_types_seen[seen_gene]} genes)')

    log.info(f'Added {added:,} features')

    gl = genelist()
    gl.load_list(promoters)
    gl.save(f'{script_path}/../../genome/{genome}_promoters.glb')

    gl_tes = genelist()
    gl_tes.load_list(newl)
    gl_tes.save(f'{script_path}/../../genome/{genome}_tes.glb') # Includes genes and TEs;

    log.info(f'Genome annotation has {len(gl_tes):,} features in total')

    print(gl_tes)

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



