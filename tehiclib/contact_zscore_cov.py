##=========================================
## normalize by total peak length of chip-seq data
##=========================================

# Once you get this working, convert into a script that takes a BED file, and adds it to the
# Figure 3A plot.

import os
import sys
import glob
import pickle
import numpy
from scipy.stats import zscore
#import pandas as pd
from matplotlib import pyplot as plot

class contact_z_score_cov:
    def __init__(self):
        self.data = self.load_data()

    def load_bed(self,
                 bed_filename:str,
                 hic_data_BEDPE:str
                 ):
        """
        **Purpose**
            Main entry point. Provide your BED file of genome locations here;

            This can be called multiple times to add multiple BED files;

        **Arguments**
            bed_filename (Required)
                The name of the BED file to be loaded.the BED file containing the locations of your peaks;

            hic_data_BEDPE (Required)
                The BEDPE file containing your list 
        """
        pass

    def load_data(self) -> dict:
        with open('../data/all_data.pkl', 'rb') as data_filename:
            data = pickle.load(data_filename)
        return data

    def _plot_peak_length_histo(self):
        import seaborn as sns

        ## distribution of peak length
        plt.figure(figsize=(9,4), dpi=300)
        sns.displot(data=peaklen,bins=40,alpha=0.5,kde=True)
        plt.xlabel('Total peak length')
        plt.tight_layout()
        plt.savefig('human_peaklen_hist.pdf')
        plt.close()

    def calc_contact_Z(self):
        """
        **Purpose**
        Calculate contact_Z versus peak length data, for background comparison versus your
        BED file;
        """
        self.chip_data_names = sorted(self.data['bkgds'].keys())
        self.peaklens = [self.data['peaklens'][chip] for chip in self.chip_data_names]
        ## mean of the background
        t_backgrnd = numpy.array([self.data['bkgds'][chip][10:15] for chip in self.chip_data_names])
        t_contacts = numpy.array([self.data['reals'][chip][10:15] for chip in self.chip_data_names])
        # normalization, (real-pseudo)/total length
        normed = t_contacts - t_backgrnd
        peaklens = numpy.array([self.data['peaklens'][chip] for chip in self.chip_data_names])
        normed = normed / peaklens[:, numpy.newaxis]

        ## zscore
        self.zscore = zscore(normed)
        self.zscore = numpy.mean(self.zscore, axis=1)

    def plot_contact_Z_scatter(self, filename):
        """
        **Purpose**

        """
        #if not self.done:
        #    print('No added BED file. Analysis is incomplete')

        fig = plot.figure()
        ax = fig.add_subplot(111)

        #print(len(self.peaklens), self.peaklens)
        #print(len(self.zscore), self.zscore)

        ax.scatter(self.peaklens, self.zscore, alpha=0.3, color='lightgrey')

        # plot a few landmarks:
        lands = set(['SMARCA4', 'MAFK', 'EZH2', 'KDM4A', 'CTCF', 'RAD21'])
        for n, x, y in zip(self.chip_data_names, self.peaklens, self.zscore):
            n = n.split('_')[0]
            if n in lands:
                ax.text(x, y, n, ha='center', va='center', fontsize=6)

        fig.savefig(filename)

if __name__ == '__main__':
    # tester.
    cov_z = contact_z_score_cov()
    cov_z.calc_contact_Z()
    #cov_z.load_bed('../data/all_data.')
    cov_z.plot_contact_Z_scatter('test.pdf')