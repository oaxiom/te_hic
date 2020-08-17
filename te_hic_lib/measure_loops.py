
import gzip, numpy

class measure_loops:
    def __init__(self, logger):
        self.logger = logger

    def bed_to_bed(self,
        reads,
        bed,
        window,
        outfile,
        **kargs):
        '''
        Measure the loops between a BED file

        '''

        if '.gz' in reads:
            readsin = gzip.open(reads, 'rt')
        else:
            readsin = open(reads, 'rt')

        if '.gz' in bed:
            bedin = gzip.open(bed, 'rt')
        else:
            bedin = open(bed, 'rt')

        # read all the BED peaks in, and set up the storage container
        peaks = {}
        for len_peaks, peak in enumerate(bedin):
            peak = peak.strip().split('\t')

            chrom = peak[0]

            if chrom not in peaks:
                peaks[chrom] = []

            # We bin the peak to the nearest window
            cpt = (int(peak[1]) + int(peak[2])) // 2
            bin_left = (cpt // window) * window

            loc = bin_left

            peaks[chrom].append(loc)
        bedin.close()

        self.logger.info('Found {0:,} BED peaks'.format(len_peaks))

        store = {}

        for idx, read_pair in enumerate(readsin):
            read_pair = read_pair.strip().split('\t')

            chrom_left = read_pair[0]
            chrom_rite = read_pair[3]

            cpt = (int(read_pair[1]) + int(read_pair[2])) // 2
            bin_left = (cpt // window) * window
            cpt = (int(read_pair[4]) + int(read_pair[5])) // 2
            bin_rite = (cpt // window) * window

            try:
                # Found a loop between two peaks in the BED
                if bin_left in peaks[chrom_left] and bin_rite in peaks[chrom_rite]:
                    korder = sorted([(chrom_left, bin_left), (chrom_rite, bin_rite)])
                    korder = korder[0] + korder[1]

                    if korder not in store:
                        store[korder] = 0
                    store[korder] += 1
            except KeyError:
                pass # chrom is not in peaks, but is in reads

            if (idx+1) % 1e6 == 0:
                self.logger.info('{0:,} reads processed'.format(idx+1))
                #break

        readsin.close()

        # Work out histogram;
        hist_max = 21
        all_scores = list(store.values())
        h = numpy.histogram(all_scores, range=[1,hist_max], bins=hist_max-1)
        self.logger.info('Histogram of loops:')
        tot = sum(i for i in h[0])
        perc = [i/tot*100 for i in h[0]]
        for v, b, p in zip(h[0], h[1], perc):
            if b == hist_max-1:
                self.logger.info('  {0}+= {1} ({2:.1f}%) reads'.format(int(b), v, p))
            else:
                self.logger.info('  {0} = {1} ({2:.1f}%) reads'.format(int(b), v, p))

        oh = gzip.open(outfile, 'wt')
        oh.write('\t'.join(['chrom1', 'left1', 'right1', 'chrom2', 'left1', 'right2', 'read_count']))
        for loop in store:
            oh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(loop[0], loop[1], loop[1]+window,
                loop[2], loop[3], loop[3]+window,
                store[loop]))

        self.logger.info('Saved {0}'.format(outfile))

        oh.close()

        return
