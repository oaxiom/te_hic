#!/usr/bin/env python3
'''

This one takes the output from assign_to_te.py
and meausres things likt
TE -> TE
TE -> -
- -> -

and collects a variety of info about the sort of linkes observed.

'''

import sys, os, gzip
from itertools import product
import miniglbase3
import common

class quantify:
    def __init__(self, project_name):
        self.project_name = project_name

    def bind_genome(self, genome_glb):
        self.genome = glbase3.glload(genelist_glb_filename)

    def measure_te_anchors(self):
        '''
        **Purpose**
            Make a crude measure of the
            TE <-> TE
            TE <-> -
            -  <-> -

            possible arrangements

        '''
        te = {'TE <-> TE': 0,
            'TE <-> -': 0,
            '-  <-> -': 0}
        res_te_te = {}
        res_te_nn = {}

        # The format of the TE file is:
        # read1.chrom read1.left read1.right read1.labels read1.type read2.chrom read2.left read2.right read2.labels read2.type

        print("Measures anchors...")
        total = 0
        oh = gzip.open(self.filename, 'rt')
        for line in oh:
            r = line.strip().split('\t')
            # measure TE anchors
            if 'TE' in r[4] and 'TE' in r[9]:
                te['TE <-> TE'] += 1
            elif 'TE' in r[4] or 'TE' in r[9]:
                te['TE <-> -'] += 1
            else:
                te['-  <-> -'] += 1
            total += 1

            if total % 1000000 == 0:
                print('Processed: {:,}'.format(total))
                #break

            # Measure TEs in detail
            if 'TE' in r[4] and 'TE' in r[9]:
                # possible to have more than one TE:
                tel = [i.strip() for i in r[3].split(',') if ':' in i] # can also hoover up some genes, so use ':' to discriminate TEs
                ter = [i.strip() for i in r[8].split(',') if ':' in i]
                combs = product(tel, ter)
                combs = [tuple(sorted(i)) for i in combs] # sort to make it unidirectional

                combs = set(combs)
                for c in combs:
                    if c not in res_te_te:
                        res_te_te[c] = 0
                    res_te_te[c] += 1

            elif 'TE' in r[4] or 'TE' in r[9]:
                if 'TE' in r[4]:
                    TE = [i.strip() for i in r[3].split(',') if ':' in i]
                elif 'TE' in r[9]:
                    TE = [i.strip() for i in r[8].split(',') if ':' in i]
                for t in TE:
                    if ':' not in t:
                        continue
                    if t not in res_te_nn:
                        res_te_nn[t] = 0
                    res_te_nn[t] += 1

        oh.close()

        print('\nmeasure_te_anchors():')
        print('  TE <-> TE : {:,} ({:.2%})'.format(te['TE <-> TE'], te['TE <-> TE']/total))
        print('  TE <-> -- : {:,} ({:.2%})'.format(te['TE <-> -'], te['TE <-> -']/total))
        print('  -- <-> -- : {:,} ({:.2%})'.format(te['-  <-> -'], te['-  <-> -']/total))
        print()

        oh = open('%s_crude_measures.txt' % self.project_name, 'w')
        oh.write('TE <-> TE : {:,} ({:.5%})\n'.format(te['TE <-> TE'], te['TE <-> TE']/total))
        oh.write('TE <-> -- : {:,} ({:.5%})\n'.format(te['TE <-> -'], te['TE <-> -']/total))
        oh.write('-- <-> -- : {:,} ({:.5%})\n'.format(te['-  <-> -'], te['-  <-> -']/total))
        oh.close()

        oh_te_te = open('%s_te-te_anchor_frequencies.tsv' % self.project_name, 'w')
        oh_te_te.write('%s\n' % '\t'.join(['TE1', 'TE2', 'RPM', 'RPM per kbp TE', 'TE1_genome_freq', 'TE2_genome_freq']))
        for k in sorted(list(res_te_te)):
            te1 = self.genome._findDataByKeyLazy('name', k[0])
            te2 = self.genome._findDataByKeyLazy('name', k[1])
            rpm = res_te_te[k]/total*1e6
            joint_kb_size = te1['genome_count'] + te2['genome_count']
            rpmpkbte = (rpm / joint_kb_size)*1e3
            line = {'te1': k[0], 'te2': k[1],
                'rpm': rpm,
                'rpmpkbte': rpmpkbte,
                'te1_genome_freq': te1['genome_percent'] / 100.0, # Convert back to fraction;
                'te2_genome_freq': te2['genome_percent'] / 100.0,
                #'enrichment': #!?!?!
                }
            oh_te_te.write('{i[te1]}\t{i[te2]}\t{i[rpm]}\t{i[rpmpkbte]}\t{i[te1_genome_freq]}\t{i[te2_genome_freq]}\n'.format(i=line))
            #print('{i[te1]}\t{i[te2]}\t{rpm}\t{rpmkbte}\t{te1_genome_freq}\t{te2_genome_freq}\n'.format(i=line))
        oh_te_te.close()

        te_nn = glbase3.genelist()
        te_nn.load_list([{'name': k, 'count': res_te_nn[k]} for k in res_te_nn])
        te_nn = te_nn.map(genelist=self.genome, key='name')
        for te in te_nn:
            te['RPM'] = (res_te_nn[te['name']]/total) * 1e6
            te['RPM per kbp of TE'] = (te['RPM'] / te['genome_count']) * 1e3
            #te['enrichment'] =
        te_nn._optimiseData()
        te_nn.sort('name')
        te_nn.saveTSV('%s_te-nn_anchor_frequencies.tsv' % self.project_name, key_order=['name', 'count', 'genome_count', 'genome_percent', 'RPM'])

        return

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('\nNot enough arguments')
        print('quantify_links.py in.tsv species project_name')
        common.print_species()
        print()
        sys.exit()

    species = sys.argv[2]
    if not common.check_species(species):
        sys.exit()

    script_path = os.path.dirname(os.path.realpath(__file__))

    q = quantify(sys.argv[3])
    q.bind_genome(os.path.join(script_path, 'genome/%s_te_genome_freqs.glb' % species))
    q.load_tsv(sys.argv[1])
    q.measure_te_anchors()


