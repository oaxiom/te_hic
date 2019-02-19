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

'''

import sys, os, math, numpy, shutil
#import glbase3
import common

class build_matrices:
    def __init__(self, script_path, species, infilename, resolution):
        '''
        Load the species and set up the bins
        '''
        # get the samplename:
        self.sample_name = os.path.split(infilename)[1].replace('.te.annot.tsv', '').replace('.tsv', '') # second is in case the user is messing with the pattern

        self.res = resolution

        self.all = {} # 1
        self.tete = {} # 2
        self.tenn = {} # 3
        self.nnnn = {} # 4

        # get the genome_sizes:
        self.chrom_sizes = {}
        self.num_bins = {}
        self.chrom_bin_offsets = {} # the top and bottom bins for this chrom

        oh = open(os.path.join(script_path, 'genome/%s.chromSizes.clean' % species), 'r')
        for line in oh:
            line = line.strip().split('\t')
            chrom = line[0]
            num_bins = math.ceil(int(line[1]) / self.res)

            self.chrom_sizes[chrom] = int(line[1])
            self.num_bins[chrom] = num_bins

        # work out the min_bin, max_bin
        cumm_bin = 0
        for chrom in sorted(self.num_bins):
            self.chrom_bin_offsets[chrom] = (cumm_bin, cumm_bin+self.num_bins[chrom])
            cumm_bin += self.num_bins[chrom]+1
        oh.close()
        return

    def build_matrices(self, infilename):
        '''

        Build the raw matrices in the style of hicpro

        TODO: Implement ICE normalisation

        And also output the hiccys?

        '''
        print('Building in-memory matrices')

        done = 0
        oh = open(infilename, 'r')
        for line in oh:
            line = line.strip().split('\t')

            # The format of the TE file is:
            # read1.chrom read1.left read1.right read1.labels read1.type read2.chrom read2.left read2.right read2.labels read2.type

            read1_chrom = line[0]
            read2_chrom = line[5]
            read1_mid = (int(line[1]) + int(line[2])) // 2
            read2_mid = (int(line[6]) + int(line[7])) // 2

            read1_bin = (read1_mid // self.res) + self.chrom_bin_offsets[read1_chrom][0]
            read2_bin = (read2_mid // self.res) + self.chrom_bin_offsets[read2_chrom][0]

            bin_pair = tuple(sorted([read1_bin, read2_bin]))

            # All:
            if bin_pair not in self.tete:
                self.all[bin_pair] = 0
            self.all[bin_pair] += 1

            if 'TE' in line[4] and 'TE' in line[9]:
                # TE <=> TE
                if bin_pair not in self.tete:
                    self.tete[bin_pair] = 0
                self.tete[bin_pair] += 1
            elif 'TE' in line[4] or 'TE' in line[9]:
                # TE <=> non-TE
                if bin_pair not in self.tenn:
                    self.tenn[bin_pair] = 0
                self.tenn[bin_pair] += 1
            else:
                # non-TE <=> non-TE
                if bin_pair not in self.nnnn:
                    self.nnnn[bin_pair] = 0
                self.nnnn[bin_pair] += 1

            done += 1
            if done % 1000000 == 0:
                print('Processed: {:,}'.format(done))
                #if done >20000000: break

        oh.close()

        return

    def save_matrices(self, out_path):
        '''
        **Purpose**
            Save the matrices into out_path/sample/resolution/
        '''
        if not os.path.isdir(out_path):
            os.mkdir(out_path) # don't delete otherwise this will be unfirendly to others working here

        if not os.path.isdir(os.path.join(out_path, self.sample_name)):
            os.mkdir(os.path.join(out_path, self.sample_name))

        if not os.path.isdir(os.path.join(out_path, self.sample_name, str(self.res))):
            os.mkdir(os.path.join(out_path, self.sample_name, str(self.res)))

        # save the BED file describing the binIDs
        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s_abs.bed' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for chrom in self.chrom_bin_offsets:
            #print(chrom)
            for localbinid, binid in enumerate(range(self.chrom_bin_offsets[chrom][0], self.chrom_bin_offsets[chrom][1])):
                l = localbinid * self.res
                r = (localbinid * self.res) + self.res
                oh.write('%s\t%s\t%s\t%s\n' % (chrom, l, r, binid+1)) # The +1 is to mimic HiCpro, remember to also +1 below!!
        oh.close()
        print('Saved BED bins: "%s"' % filename)

        # Save the matrices:
        # matrices are sparse:
        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.all.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.all):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.all[bins]))
        oh.close()
        print('Saved All matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.tete.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.tete):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.tete[bins]))
        oh.close()
        print('Saved TE <=> TE matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.tenn.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.tenn):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.tenn[bins]))
        oh.close()
        print('Saved TE <=> non-TE matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.nnnn.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.nnnn):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.nnnn[bins]))
        oh.close()
        print('Saved non-TE <=> non-TE matrix: "%s"' % filename)


if __name__ == '__main__':
    if len(sys.argv) != 5:
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

    mat = build_matrices(script_path, species, sys.argv[3], int(sys.argv[2]))
    mat.build_matrices(sys.argv[3])
    mat.save_matrices(sys.argv[4])



