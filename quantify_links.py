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
#import glbase3

class quantify:
    def __init__(self, project_name):
        self.project_name = project_name

    def load_tsv(self, filename):
        '''
        load the output from assign_to_te.py
        '''

        self.reads = []
        oh = open(filename, 'r')
        for line in oh:
            line = line.strip().split('\t')
            self.reads.append(line)
        oh.close()

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

        total = 0
        for r in self.reads:
            if 'TE' in r[4] and TE in r[8]:
                te['TE <-> TE'] += 1
            elif 'TE' in r[4] or TE in r[8]:
                te['TE <-> -'] += 1
            else:
                te['-  <-> -'] += 1
            total += 1

        print('\nmeasure_te_anchors():')
        print('  TE <-> TE : {:,} ({:.2%})'.format(te['TE <-> TE'], te['TE <-> TE']/total))
        print('  TE <-> - : {:,} ({:.2%})'.format(te['TE <-> -'], te['TE <-> -']/total))
        print('  -  <-> - : {:,} ({:.2%})'.format(te['-  <-> -'], te['-  <-> -']/total))
        print()

    def measure_te_freqs(self):
        '''
        **Purpose**
            Measure the types of TE frequencies, between pairs of TE -> TE and TE -> -
        '''
        res_te_te = {}
        res_te_nn = {}

        for idx, r in enumerate(self.reads):
            if 'TE' in r[4] and TE in r[8]:
                # possible to have more than one TE:
                tel = [i.strip() for i in r[4].split(',')]
                ter = [i.strip() for i in r[8].split(',')]
                combs = product(tel, ter)

                combs = set(combs)
                
            elif 'TE' in r[4] or TE in r[8]:
                pass

            if i > 10:
                break

            if done % 1000 == 0:
                print('Processed: {:,}'.format(done)) 
                break

        oh_te_te = open('%s_te-nn_anchor_frequencies.tsv' % self.project_name)
        oh_te_te.close()

        oh_te_nn = open('%s_te-te_anchor_frequencies.tsv' % self.project_name)
        oh_te_nn.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('\nNot enough arguments')
        print('quantify_links.py in.tsv project_name')
        print()
        sys.exit()
    
    q = quantify()
    q.load_tsv(sys.argv[1])
    q.measure_te_anchors()
    q.measure_te_freqs('%s_te_anchor_frequencies.tsv' % sys.argv[2])


