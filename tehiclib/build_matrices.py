'''

Part of the te_hic

Build four (raw) matrices:
1. All
2. TE<=>TE only
3. TE<=>non-TE only
4. non-TE <=> non-TE only

TODO:
5. TE<=>non-TE and TE<=>TE only?

'''

import sys, os, math, numpy, shutil, gzip
from . import common

class build_matrices:
    def __init__(self, genome_sizes, resolution, logger):
        '''
        Load the species and set up the bins
        '''

        self.log = logger

        # get the samplename:
        self.res = resolution

        self.all = {} # 1
        self.tete = {} # 2
        self.tenn = {} # 3
        self.nnnn = {} # 4

        # get the genome_sizes:
        self.chrom_sizes = {}
        self.num_bins = {}
        self.chrom_bin_offsets = {} # the top and bottom bins for this chrom

        for chrom in genome_sizes:
            num_bins = math.ceil(genome_sizes[chrom]) / self.res)
            self.chrom_sizes[chrom] = int(line[1])
            self.num_bins[chrom] = num_bins

        # work out the min_bin, max_bin
        cumm_bin = 0
        for chrom in sorted(self.num_bins):
            self.chrom_bin_offsets[chrom] = (cumm_bin, cumm_bin+self.num_bins[chrom])
            cumm_bin += self.num_bins[chrom]+1
        oh.close()
        return

    def build_matrices(self, mapped_pairs):
        '''

        Build the raw matrices in the style of hicpro

        TODO: Implement ICE normalisation

        And also output the hiccys?

        '''
        self.log.info('Building in-memory matrices')


        for done, pair in enumerate(mapped_pairs):
            print(pair)
            1/0

            line = line.strip().split('\t')

            read1_chrom = line[0]
            read2_chrom = line[5]
            read1_mid = (int(line[1]) + int(line[2])) // 2
            read2_mid = (int(line[6]) + int(line[7])) // 2

            read1_bin = (read1_mid // self.res) + self.chrom_bin_offsets[read1_chrom][0]
            read2_bin = (read2_mid // self.res) + self.chrom_bin_offsets[read2_chrom][0]

            bin_pair = tuple(sorted([read1_bin, read2_bin]))

            if read1_bin == read2_bin: # ignore all mid-line self-bins
                continue

            # All:
            if bin_pair not in self.all:
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
                print(f'Processed: {done:,}')

        oh.close()

        return

    def save_matrices(self, out_path):
        '''
        **Purpose**
            Save the matrices into out_path/sample/resolution/
        '''
        if not os.path.isdir(out_path):
            os.mkdir(out_path) # don't delete otherwise this will be unfriendly to others working here

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
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.all[bins])) # The +1 is to mimic HiCpro!
        oh.close()
        print('Saved All matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.tete.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.tete):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.tete[bins]))# The +1 is to mimic HiCpro!
        oh.close()
        print('Saved TE <=> TE matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.tenn.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.tenn):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.tenn[bins])) # The +1 is to mimic HiCpro!
        oh.close()
        print('Saved TE <=> non-TE matrix: "%s"' % filename)

        filename = os.path.join(out_path, self.sample_name, str(self.res), '%s_%s.nnnn.raw.matrix' % (self.sample_name, self.res))
        oh = open(filename, 'w')
        for bins in sorted(self.nnnn):
            oh.write('%s\t%s\t%s\n' % (bins[0]+1, bins[1]+1, self.nnnn[bins])) # The +1 is to mimic HiCpro!
        oh.close()
        print('Saved non-TE <=> non-TE matrix: "%s"' % filename)

