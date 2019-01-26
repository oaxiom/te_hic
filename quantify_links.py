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
        self.filename = filename

        return
        # Below is deprecated. It's fast, but when the reads gets above ~100 million it tends to kill
        # the computer even with >64G RAM
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
        res_te_te = {}
        res_te_nn = {}

        print("Measures anchors...")
        total = 0
        oh = open(self.filename, 'r')
        for line in oh:
            r = line.strip().split('\t')
            # measure TE anchors
            if 'TE' in r[4] and 'TE' in r[8]:
                te['TE <-> TE'] += 1
            elif 'TE' in r[4] or 'TE' in r[8]:
                te['TE <-> -'] += 1
            else:
                te['-  <-> -'] += 1
            total += 1
            if total % 1000000 == 0:
                print('Processed: {:,}'.format(total))
                #break

            # MEasure TEs in detail
            if 'TE' in r[4] and 'TE' in r[8]:
                # possible to have more than one TE:
                tel = [i.strip() for i in r[3].split(',') if ':' in i] # can also hoover up some genes, so use ':' to discriminate TEs
                ter = [i.strip() for i in r[7].split(',') if ':' in i]
                combs = product(tel, ter)
                combs = [tuple(sorted(i)) for i in combs] # sort to make it unidirectional                

                combs = set(combs)
                for c in combs:
                    if c not in res_te_te:
                        res_te_te[c] = 0
                    res_te_te[c] += 1

            elif 'TE' in r[4] or 'TE' in r[8]:
                if 'TE' in r[4]:
                    TE = [i.strip() for i in r[3].split(',') if ':' in i]
                elif 'TE' in r[8]:
                    TE = [i.strip() for i in r[7].split(',') if ':' in i]
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
        oh_te_te.write('%s\n' % '\t'.join(['TE1', 'TE2', 'count', '%']))
        for k in sorted(list(res_te_te)):
            oh_te_te.write('%s\t%s\t%s\t%.6f\n' % (k[0], k[1], res_te_te[k], res_te_te[k]/total*100.0))
        oh_te_te.close()

        oh_te_nn = open('%s_te-nn_anchor_frequencies.tsv' % self.project_name, 'w')
        oh_te_nn.write('%s\n' % '\t'.join(['TE1', 'count', '%']))
        for k in sorted(list(res_te_nn)):
            oh_te_nn.write('%s\t%s\t%.6f\n' % (k, res_te_nn[k], res_te_nn[k]/total))
        oh_te_nn.close()

        return

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('\nNot enough arguments')
        print('quantify_links.py in.tsv project_name')
        print()
        sys.exit()
    
    q = quantify(sys.argv[2])
    q.load_tsv(sys.argv[1])
    q.measure_te_anchors()


