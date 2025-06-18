##=========================================
## normalize by total peak length of chip-seq data
##=========================================

# Once you get this working, convert into a script that takes a BED file, and adds it to the
# Figure 3A plot.

import os
import sys
import glob
import pickle
import gzip
import numpy
import random
from scipy.stats import zscore
#import pandas as pd
from matplotlib import pyplot as plot

class contact_z_score_cov:
    def __init__(self, logger):
        self.data = self.load_data()
        self.logger = logger

    def load_data(self) -> dict:
        data_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../data/all_data.pkl')

        with open(data_filename, 'rb') as data_file:
            data = pickle.load(data_file)
        return data

    def insert_contacts(self,
                        bed_contact_list,
                        number_of_peaks,
                        peaklen,
                        label,
                        random_background,
                        GC=False,
                        ):
        """
        **Purpose**
            insert your custom peak list here;

            This can be called multiple times to add multiple BED files;
        """
        # convert the dict to list;
        real = [bed_contact_list[i] for i in sorted(bed_contact_list.keys())]
        rand = [random_background[i] for i in sorted(random_background.keys())]

        # simulate background...
        self.data['reals'][f'{label} (Insert)'] = real
        self.data['peaklens'][f'{label} (Insert)'] = peaklen
        if GC:
            self.data['gc_bkgds'][f'{label} (Insert)'] = rand
        else:
            self.data['bkgds'][f'{label} (Insert)'] = rand

        self.GC = GC

    def calc_contact_Z(self):
        """
        **Purpose**
        Calculate contact_Z versus peak length data, for background comparison versus your
        BED file;
        """
        self.chip_data_names = sorted(self.data['reals'].keys())
        self.peaklens = [self.data['peaklens'][chip] for chip in self.chip_data_names]

        ## mean of the background
        if self.GC:
            t_backgrnd = numpy.array([self.data['gc_bkgds'][chip][10:18] for chip in self.chip_data_names])
        else:
            t_backgrnd = numpy.array([self.data['bkgds'][chip][10:18] for chip in self.chip_data_names])

        t_contacts = numpy.array([self.data['reals'][chip][10:18] for chip in self.chip_data_names])
        #print(t_backgrnd, t_contacts)
        # normalization, (real-pseudo)/total length
        normed = t_contacts - t_backgrnd
        peaklens = numpy.array([self.data['peaklens'][chip] for chip in self.chip_data_names])
        normed = normed / peaklens[:, numpy.newaxis]

        ## zscore
        self.zscore = zscore(normed)
        self.zscore = numpy.mean(self.zscore, axis=1)

    def __quartiles(self, data):
        # Calculate quartiles
        Q1 = numpy.quantile(data, 0.25)
        Q3 = numpy.quantile(data, 0.75)

        return Q1, Q3

    def plot_contact_Z_scatter(self, filename):
        """
        **Purpose**

        """
        fig = plot.figure()
        ax = fig.add_subplot(111)

        # quartiles
        Q1, Q3 = self.__quartiles(self.zscore)
        self.logger.info(f"First Quartile ({Q1:.2f}):")
        self.logger.info(f"Third Quartile ({Q3:.2f}):")

        formers = []
        breakers = []

        spot_cols = []
        for n, x, y in zip(self.chip_data_names, self.peaklens, self.zscore):
            if '(Insert)' in n:
                spot_cols.append('tab:red')
                self.logger.info(f'The predicted contact Z-score for {n} is {y:.2f}')
                ax.set_title(f'The predicted contact Z-score for {n.replace("(Insert)", "")} is {y:.2f}')
            elif y >= Q3:
                spot_cols.append('tab:green')
                formers.append((n, x))
            elif y <= Q1:
                spot_cols.append('tab:blue')
                breakers.append((n, x))
            else:
                spot_cols.append('lightgrey')

        ax.scatter(self.peaklens, self.zscore, alpha=0.3, color=spot_cols)

        # plot a few landmarks:
        lands = set(['SMARCA4', 'KDM4A', 'CTCF', 'RAD21'])
        for n, x, y in zip(self.chip_data_names, self.peaklens, self.zscore):
            split_name = n.split('_')[0]
            if split_name in lands:
                ax.text(x, y, split_name, ha='center', va='center', fontsize=6)
            elif '(Insert)' in n:
                ax.text(x, y, n, ha='center', va='center', fontsize=6)

        ax.axhline(0, lw=0.5, color='grey')
        ax.axhline(Q3, ls=':', lw=0.5, color='grey')
        ax.axhline(Q1, ls=':', lw=0.5, color='grey')

        ax.set_ylabel('Contact Z-score')
        ax.set_xlabel('Number of bp in peaks')

        fig.savefig(filename)

        return formers, breakers

    def __load_bed(self, filename):
        # Not the same as measure_contacs.__load_bed
        if '.gz' in filename:
            bedin = gzip.open(filename, 'rt')
        else:
            bedin = open(filename, 'rt')

        # read all the BED peaks in, and set up the storage container
        peaklen = 0 # number of base pairs occupied by the peaks;
        peaks = {}
        for len_peaks, peak in enumerate(bedin):
            peak = peak.strip().split('\t')

            chrom = peak[0]

            if chrom not in peaks:
                peaks[chrom] = []

            l = int(peak[1])
            r = int(peak[2])
            peaklen += r - l

            peaks[chrom].append((l, r))
        bedin.close()

        return peaks, peaklen, len_peaks

    def __load_bed_GC(self, filename):
        # Not the same as measure_contacs.__load_bed
        if '.gz' in filename:
            bedin = gzip.open(filename, 'rt')
        else:
            bedin = open(filename, 'rt')

        # read all the BED peaks in, and set up the storage container
        peaklen = 0 # number of base pairs occupied by the peaks;
        peaks = {}
        for len_peaks, peak in enumerate(bedin):
            peak = peak.strip().split('\t')

            chrom = peak[0]

            if chrom not in peaks:
                peaks[chrom] = []

            l = int(peak[1])
            r = int(peak[2])
            gc = int(peak[4])
            peaklen += r - l

            peaks[chrom].append((l, r, gc))
        bedin.close()

        return peaks, peaklen, len_peaks

    def __generate_matched_random_GC(self, bed_file) -> dict:
        """
        **Emulate bedtools shuf, but without a genome file;
        and match GC content
        """
        peaks, peak_len_in_bp, len_peaks = self.__load_bed_GC(bed_file)

        rand_peaks = {}
        for chrom in peaks:
            rand_peaks[chrom] = []
            #chrom = chrom.lstrip('chr')
            #print(chrom)

            for peak in peaks[chrom]:
                # get a peak from the background randon;
                gc = peak[2]

                # There are too few loci to get full coverage at these extreme GC contents, so bracket them.
                if gc < 20: gc = 20
                if gc > 70: gc = 70

                rand_peak = random.choice(self.data['randoms_gc'][gc][chrom])
                # resize to same size
                psz = peak[1] - peak[0]

                rand_peak = (rand_peak, rand_peak + psz)
                rand_peaks[chrom].append(rand_peak)

        # check the new_peak_len_in_bp
        new_peak_len_in_bp = 0
        new_len_peaks = 0
        for chrom in rand_peaks:
            for peak in rand_peaks[chrom]:
                new_peak_len_in_bp += peak[1] - peak[0]

            new_len_peaks += len(rand_peaks[chrom])

        self.logger.info('Shuffled BED, sanity check:')
        self.logger.info(f'Number of bp in peaks: {peak_len_in_bp} = {new_peak_len_in_bp}')
        self.logger.info(f'Number of peaks: {len_peaks} = {new_len_peaks-1}')

        return dict(peaks=rand_peaks, peak_len_in_bp=new_peak_len_in_bp, len_peaks=new_len_peaks)

    def generate_matched_random(self, bed_file, GC=False) -> dict:
        """
        **Emulate bedtools shuf, but without a genome file;
        """
        if GC:
            return self.__generate_matched_random_GC(bed_file)

        peaks, peak_len_in_bp, len_peaks = self.__load_bed(bed_file)

        rand_peaks = {}
        for chrom in peaks:
            rand_peaks[chrom] = []

            for peak in peaks[chrom]:
                # get a peak from the background randon;
                try:
                    rand_peak = random.choice(self.data['randoms'][chrom])
                except KeyError: # Bad chrom;
                    rand_peak = random.choice(self.data['randoms']['chr1'])
                # resize to same size
                psz = peak[1] - peak[0]
                rand_peak = (rand_peak, rand_peak + psz)
                rand_peaks[chrom].append(rand_peak)

        # check the new_peak_len_in_bp
        new_peak_len_in_bp = 0
        new_len_peaks = 0
        for chrom in rand_peaks:
            for peak in rand_peaks[chrom]:
                new_peak_len_in_bp += peak[1] - peak[0]

            new_len_peaks += len(rand_peaks[chrom])

        self.logger.info('Shuffled BED, sanity check:')
        self.logger.info(f'Number of bp in peaks: {peak_len_in_bp} = {new_peak_len_in_bp}')
        self.logger.info(f'Number of peaks: {len_peaks} = {new_len_peaks}')

        return dict(peaks=rand_peaks, peak_len_in_bp=new_peak_len_in_bp, len_peaks=new_len_peaks)

if __name__ == '__main__':
    # tester.
    cov_z = contact_z_score_cov()
    cov_z.calc_contact_Z()
    #cov_z.load_bed('../data/all_data.')
    cov_z.plot_contact_Z_scatter('test.pdf')