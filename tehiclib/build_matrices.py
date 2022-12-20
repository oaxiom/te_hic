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
        self.res = resolution * 1000

        self.all = {} # 1
        self.tete = {} # 2
        self.tenn = {} # 3
        self.nnnn = {} # 4

        # get the genome_sizes:
        self.chrom_sizes = {}
        self.num_bins = {}
        self.chrom_bin_offsets = {} # the top and bottom bins for this chrom

        for chrom in genome_sizes:
            num_bins = math.ceil(genome_sizes[chrom] / self.res)
            self.num_bins[chrom] = num_bins

        # work out the min_bin, max_bin
        cumm_bin = 0
        for chrom in sorted(self.num_bins):
            self.chrom_bin_offsets[chrom] = (cumm_bin, cumm_bin+self.num_bins[chrom])
            cumm_bin += self.num_bins[chrom]+1

        return

    def build_matrices(self, mapped_pairs):
        '''

        Build the raw matrices in the style of hicpro

        '''
        self.log.info(f'Building in-memory matrices for resolution {self.res} bp')

        # The format of mapped_pairs is:
        # (
        #   ('8', 88996379, '8', 89417593, '+', '-'),
        #   {'Tigger19a:TcMar-Tigger:DNA'}, {'TE'},
        #   {'MER33:hAT-Charlie:DNA'}, {'TE'}
        # )

        for done, pair in enumerate(mapped_pairs):
            read1_chrom = f'chr{pair[0][0]}'
            read2_chrom = f'chr{pair[0][2]}'
            read1_mid = pair[0][1]
            read2_mid = pair[0][3]

            read1_bin = (read1_mid // self.res) + self.chrom_bin_offsets[read1_chrom][0]
            read2_bin = (read2_mid // self.res) + self.chrom_bin_offsets[read2_chrom][0]

            bin_pair = tuple(sorted([read1_bin, read2_bin]))

            if read1_bin == read2_bin: # ignore all mid-line self-bins
                continue

            # All:
            if bin_pair not in self.all:
                self.all[bin_pair] = 0
            self.all[bin_pair] += 1

            if 'TE' in pair[2] and 'TE' in pair[4]:
                # TE <=> TE
                if bin_pair not in self.tete:
                    self.tete[bin_pair] = 0
                self.tete[bin_pair] += 1

            elif 'TE' in pair[2] or 'TE' in pair[4]:
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
            if done % 1e6 == 0:
                self.log.info(f'Processed: {done:,}')

        return

    def save_matrices(self, label):
        '''
        **Purpose**
            Save the matrices into out_path/sample/resolution/
        '''
        if not os.path.isdir(f'matrices_{label}'):
            os.mkdir(f'matrices_{label}') # don't delete otherwise this will be unfriendly to others working here

        if not os.path.isdir(os.path.join(f'matrices_{label}', str(self.res))):
            os.mkdir(os.path.join(f'matrices_{label}', str(self.res)))

        # save the BED file describing the binIDs
        filename = os.path.join(f'matrices_{label}', str(self.res), f'{label}_{self.res}_abs.bed')
        oh = open(filename, 'w')
        for chrom in self.chrom_bin_offsets:
            #print(chrom)
            for localbinid, binid in enumerate(range(self.chrom_bin_offsets[chrom][0], self.chrom_bin_offsets[chrom][1])):
                l = localbinid * self.res
                r = (localbinid * self.res) + self.res
                oh.write(f'{chrom}\t{l}\t{r}\t{binid+1}\n') # The +1 is to mimic HiCpro, remember to also +1 below!!
        oh.close()
        self.log.info(f'Saved BED bins: "{filename}"')

        # Save the matrices:
        # matrices are sparse:
        filename = os.path.join(f'matrices_{label}',str(self.res), f'{label}_{self.res}.all.raw.matrix')
        oh = open(filename, 'w')
        for bins in sorted(self.all):
            oh.write(f'{bins[0]+1}\t{bins[1]+1}\t{self.all[bins]}\n') # The +1 is to mimic HiCpro!
        oh.close()
        self.log.info(f'Saved All matrix: "{filename}"')

        filename = os.path.join(f'matrices_{label}', str(self.res), f'{label}_{self.res}.tete.raw.matrix')
        oh = open(filename, 'w')
        for bins in sorted(self.tete):
            oh.write(f'{bins[0]+1}\t{bins[1]+1}\t{self.tete[bins]}\n')# The +1 is to mimic HiCpro!
        oh.close()
        self.log.info(f'Saved TE <=> TE matrix: "{filename}"')

        filename = os.path.join(f'matrices_{label}', str(self.res), f'{label}_{self.res}.tenn.raw.matrix')
        oh = open(filename, 'w')
        for bins in sorted(self.tenn):
            oh.write(f'{bins[0]+1}\t{bins[1]+1}\t{self.tenn[bins]}\n') # The +1 is to mimic HiCpro!
        oh.close()
        self.log.info(f'Saved TE <=> non-TE matrix: "{filename}"')

        filename = os.path.join(f'matrices_{label}', str(self.res), f'{label}_{self.res}.nnnn.raw.matrix')
        oh = open(filename, 'w')
        for bins in sorted(self.nnnn):
            oh.write(f'{bins[0]+1}\t{bins[1]+1}\t{self.nnnn[bins]}\n') # The +1 is to mimic HiCpro!
        oh.close()
        self.log.info(f'Saved non-TE <=> non-TE matrix: "{filename}"')

