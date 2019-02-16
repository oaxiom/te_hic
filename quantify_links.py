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
import glbase3

class quantify:
    def __init__(self, project_name):
        self.project_name = project_name

    def bind_genome(self, genelist_glb_filename):
        self.genome = glbase3.glload(genelist_glb_filename)
        print('Loaded %s' % genelist_glb_filename)

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
            oh_te_te.write('%s\t%s\t%s\t%s\n' % (k[0], k[1], res_te_te[k], res_te_te[k]/total*100.0))
        oh_te_te.close()

        # How are you supposed to work this out?
        #te_te = glbase3.genelist()
        #te_te.load_list([{'name': str(k), 'count': res_te_nn[k]} for k in res_te_nn])
        #te_te.map((genelist=self.genome, key='name')
        #for te_pair in te_te:
        #    
        #te_te._optimiseData()
        #te_te.sort('name')
        #te_te.saveTSV('%s_te-te_anchor_frequencies.tsv' % self.project_name)    

        te_nn = glbase3.genelist()
        te_nn.load_list([{'name': k, 'count': res_te_nn[k]} for k in res_te_nn])
        te_nn = te_nn.map(genelist=self.genome, key='name')
        for te in te_nn:
            te['percent'] = (res_te_nn[te['name']]/total) * 100.0
            te['RPM'] = (res_te_nn[te['name']]/total) * 1e6
            te['RPM per kbp of TE'] = te['RPM'] / te['genome_count'] * 1e3
            #te['enrichment'] = 
        te_nn._optimiseData() 
        te_nn.sort('name')
        te_nn.saveTSV('%s_te-nn_anchor_frequencies.tsv' % self.project_name, key_order=['name', 'count', 'genome_count', 'percent', 'genome_percent', 'RPM'])
        
        #oh_te_nn = open('%s_te-nn_anchor_frequencies.tsv' % self.project_name, 'w')
        #oh_te_nn.write('%s\n' % '\t'.join(['TE1', 'count', '%']))
        #for k in sorted(list(res_te_nn)):
        #    oh_te_nn.write('%s\t%s\t%.6f\n' % (k, res_te_nn[k], res_te_nn[k]/total))
        #oh_te_nn.close()

        return

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('\nNot enough arguments')
        print('quantify_links.py in.tsv species project_name')
        print('  Valid Species codes are:')
        print('    hg38 - human')
        print('    mm10 - mouse')
        print()
        sys.exit()
    
    species = sys.argv[2]
    if species not in ('mm10', 'hg38'):
        print('Species "%s" not found' % species)
        print('Valid Species codes are:')
        print('    hg38 - human')
        print('    mm10 - mouse')
        print()
        sys.exit()

    script_path = os.path.dirname(os.path.realpath(__file__))

    q = quantify(sys.argv[3])
    q.bind_genome(os.path.join(script_path, 'genome/%s_te_genome_freqs.glb' % species))
    q.load_tsv(sys.argv[1])
    q.measure_te_anchors()


