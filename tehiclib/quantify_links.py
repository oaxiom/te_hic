#!/usr/bin/env python3
'''

This one takes the output from assign_to_te.py
and meausres things likt
TE -> TE
TE -> -
- -> -

and collects a variety of info about the sort of linkes observed.

'''

import sys, os
from itertools import product
from . import miniglbase3
from . import common

class quantify:
    def __init__(self, project_name, logger=None):
        assert logger, 'Need to provide a log'
        self.project_name = project_name
        self.log = logger

    def bind_te_freqs(self, te_genome_freqs_glb):
        self.te_freqs = te_genome_freqs_glb

    def measure_te_anchors(self, mapped_pairs_temp_file):
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

        # The format of mapped_pairs is:
        # (
        #   ('8', 88996379, '8', 89417593, '+', '-'),
        #   {'Tigger19a:TcMar-Tigger:DNA'}, {'TE'},
        #   {'MER33:hAT-Charlie:DNA'}, {'TE'}
        # )

        self.log.info("Measures anchors...")
        total = 0

        mapped_pairs = open(mapped_pairs_temp_file, 'r')

        for r in mapped_pairs:
            # (chromA, leftA, riteA, chromB, leftB, riteB, read1_feat, read1_type, read2_feat, read2_type)
            r = r.strip().split('\t')
            # measure TE anchors
            if 'TE' in r[7] and 'TE' in r[9]:
                te['TE <-> TE'] += 1
            elif 'TE' in r[7] or 'TE' in r[9]:
                te['TE <-> -'] += 1
            else:
                te['-  <-> -'] += 1
            total += 1

            if total % 1e6 == 0:
                self.log.info(f'Processed: {total:,}')
                #break

            # Measure TEs in detail
            if 'TE' in r[7] and 'TE' in r[9]:
                # possible to have more than one TE:
                tel = [i for i in r[1] if ':' in i] # can also hoover up some genes, so use ':' to discriminate TEs
                ter = [i for i in r[3] if ':' in i]
                combs = product(tel, ter)
                combs = [tuple(sorted(i)) for i in combs] # sort to make it unidirectional

                combs = set(combs)
                for c in combs:
                    if c not in res_te_te:
                        res_te_te[c] = 0
                    res_te_te[c] += 1

            elif 'TE' in r[7] or 'TE' in r[9]:
                if 'TE' in r[7]:
                    TE = [i for i in eval(r[6]) if ':' in i]
                elif 'TE' in r[9]:
                    TE = [i for i in eval(r[8]) if ':' in i]

                for t in TE:
                    if ':' not in t:
                        continue
                    if t not in res_te_nn:
                        res_te_nn[t] = 0
                    res_te_nn[t] += 1

        mapped_pairs.close()

        self.log.info('  TE <-> TE : {:,} ({:.2%})'.format(te['TE <-> TE'], te['TE <-> TE']/total))
        self.log.info('  TE <-> -- : {:,} ({:.2%})'.format(te['TE <-> -'], te['TE <-> -']/total))
        self.log.info('  -- <-> -- : {:,} ({:.2%})'.format(te['-  <-> -'], te['-  <-> -']/total))

        oh = open('stage3.%s_crude_measures.txt' % self.project_name, 'w')
        oh.write('TE <-> TE : {:,} ({:.5%})\n'.format(te['TE <-> TE'], te['TE <-> TE']/total))
        oh.write('TE <-> -- : {:,} ({:.5%})\n'.format(te['TE <-> -'], te['TE <-> -']/total))
        oh.write('-- <-> -- : {:,} ({:.5%})\n'.format(te['-  <-> -'], te['-  <-> -']/total))
        oh.close()

        oh_te_te = open(f'stage3.{self.project_name}_te-te_anchor_frequencies.tsv', 'w')
        oh_te_te.write('%s\n' % '\t'.join(['TE1', 'TE2', 'RPM', 'RPM per kbp TE', 'TE1_genome_freq', 'TE2_genome_freq']))
        for k in sorted(list(res_te_te)):
            te1 = self.te_freqs._findDataByKeyLazy('name', k[0])
            te2 = self.te_freqs._findDataByKeyLazy('name', k[1])
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
            oh_te_te.write(f'{line["te1"]}\t{line["te2"]}\t{line["rpm"]}\t{line["rpmpkbte"]}\t{line["te1_genome_freq"]}\t{line["te2_genome_freq"]}\n')

        oh_te_te.close()

        te_nn = miniglbase3.genelist()
        te_nn.load_list([{'name': k, 'count': res_te_nn[k]} for k in res_te_nn])
        te_nn = te_nn.map(genelist=self.te_freqs, key='name', silent=True)
        for te in te_nn:
            te['RPM'] = (res_te_nn[te['name']]/total) * 1e6
            te['RPM per kbp of TE'] = (te['RPM'] / te['genome_count']) * 1e3
            #te['enrichment'] =

        te_nn._optimiseData()
        te_nn.sort('name')
        te_nn.saveTSV(f'stage3.{self.project_name}_te-nn_anchor_frequencies.tsv', key_order=['name', 'count', 'genome_count', 'genome_percent', 'RPM'])

        return


