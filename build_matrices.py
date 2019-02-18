#!/usr/bin/env python3
'''

Part of the te_hic suite

Build five (raw) matrices:
1. All
2. TE<=>TE only
3. TE<=>non-TE only
4. non-TE <=> non-TE only
5. TE<=>non-TE and TE<=>TE only

'''

import sys, os, math
import glbase3
import common

def build_matrices(resolution, infilename, outpath):
    '''

    Build the raw matrices in the style of hicpro

    TODO: Implement ICE normalisation

    And also output the hiccys?

    '''
    # set up the bins



if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('\nNot enough arguments')
        print('build_matrices.py species resolution in.te.annot.tsv matrix_path')
        print('  The input file is the tsv produced by assign_to_te.py')
        print('  resolution is in base pairs (i.e. 150000)')
        print()
        print('  Valid Species codes are:')
        print('    hg38 - human')
        print('    mm10 - mouse')
        print()
        sys.exit()

    species = sys.argv[1]
    if not common.check_species(species):
        sys.exit()

    script_path = os.path.dirname(os.path.realpath(__file__))

    build_matrices(sys.argv[1], sys.argv[2], sys.argv[3])

    #mte = measureTE(sys.argv[0])
    #mte.bind_genome(os.path.join(script_path, 'genome/%s_glb_gencode_tes.glb' % species))
    #mte.load_bedpe(sys.argv[2], sys.argv[3])


