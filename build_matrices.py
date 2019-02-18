#!/usr/bin/env python3
'''

Part of the te_hic suite

Build five (raw) matrices:
1. All
2. TE<=>TE only
3. TE<=>non-TE only
4. non-TE <=> non-TE only

TODO:
5. TE<=>non-TE and TE<=>TE only
# Only supports inter-chromosome

'''

import sys, os, math, numpy
import glbase3
import common

class build_matrices:
    def __init__(self, script_path, species, resolution):
        '''
        Load the species and set up the bins
        '''
        self.res = resolution

        self.all = {} # 1
        self.tete = {} # 2
        self.tenn = {} # 3
        self.nnnn = {} # 4

        # get the genome_sizes:
        self.chrom_sizes = {}
        self.num_bins = {}

        oh = open(os.path.join(script_path, '%s.chromSizes.clean' % species), 'r')
        for line in oh:
            line = line.strip().split('\t')
            chrom = line[0]
            num_bins = math.ceil(int(line[1]) / self.res)

            self.chrom_sizes[chrom] = int(line[1])
            self.bin_sizes[chrom] = num_bins

            self.all[chrom]  = numpy.zeros(num_bins,num_bins)
            self.tete[chrom] = numpy.zeros(num_bins,num_bins)
            self.tenn[chrom] = numpy.zeros(num_bins,num_bins)
            self.nnnn[chrom] = numpy.zeros(num_bins,num_bins)

    def build_matrices(infilename):
        '''

        Build the raw matrices in the style of hicpro

        TODO: Implement ICE normalisation

        And also output the hiccys?

        '''
        oh = open(infilename, 'r')
        for line in oh:
            line = line.strip().split('\t')





if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('\nNot enough arguments')
        print('build_matrices.py species resolution in.te.annot.tsv matrix_path')
        print('  The input file is the tsv produced by assign_to_te.py')
        print('  resolution is in base pairs (i.e. 150000)')
        print()
        common.print_species()
        print()
        sys.exit()

    species = sys.argv[1]
    if not common.check_species(species):
        sys.exit()

    script_path = os.path.dirname(os.path.realpath(__file__))

    mat = build_matrices(script_path, species, int(sys.argv[2]))
    mat.build_matrices(sys.argv[2], sys.argv[3])
    # mat.save_matrices()?



