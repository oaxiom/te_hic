
import gzip, numpy

class measure_contacts:
    def __init__(self, logger):
        self.logger = logger

    def __load_bed(self, filename, window):
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

            # We bin the peak to the nearest window
            l = int(peak[1])
            r = int(peak[2])
            cpt = (l + r) // 2
            bin_left = (cpt // window) * window

            peaklen += r - l

            loc = bin_left

            peaks[chrom].append(loc)
        bedin.close()

        return peaks, peaklen, len_peaks

    def bed_to_bed(self,
        reads,
        bed,
        window,
        outfile,
        threshold: int = 1, # in reads;
        __silent:bool = False,
        **kargs):
        '''
        Measure the contacts between a BED file

        '''
        if isinstance(bed, str):
            peaks, peaklen, len_peaks = self.__load_bed(bed, window)
        elif isinstance(bed, dict):
            peaks = bed['peaks']
            peaklen = bed['peak_len_in_bp']
            len_peaks = bed['len_peaks']
        else:
            raise AssertionError('bed is not a filename or a dict')

        if not __silent: self.logger.info('Found {0:,} BED peaks'.format(len_peaks))

        if '.gz' in reads:
            readsin = gzip.open(reads, 'rt')
        else:
            readsin = open(reads, 'rt')

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
                # Found a contact between two peaks in the BED
                if bin_left in peaks[chrom_left] and bin_rite in peaks[chrom_rite]:
                    korder = sorted([(chrom_left, bin_left), (chrom_rite, bin_rite)])
                    korder = korder[0] + korder[1] # quple;

                    if korder not in store:
                        store[korder] = 0
                    store[korder] += 1
            except KeyError:
                pass # chrom is not in peaks, but is in reads

            if (idx+1) % 1e6 == 0:
                if not __silent: self.logger.info('{0:,} reads processed'.format(idx+1))

        readsin.close()

        # Work out histogram;
        hist_max = 21
        all_scores = list(store.values())
        h = numpy.histogram(all_scores, range=[1,hist_max], bins=hist_max-1)
        if not __silent: self.logger.info('Histogram of contacts:')
        tot = sum(i for i in h[0])
        perc = [i/tot*100 for i in h[0]]

        for v, b, p in zip(h[0], h[1], perc):
            if int(b) == hist_max-1:
                if int(b) == 1:
                    if not __silent: self.logger.info('  {1} ({2:.1f}%) contacts have {0}+ read'.format(int(b), v, p))
                else:
                    if not __silent: self.logger.info('  {1} ({2:.1f}%) contacts have {0}+ reads'.format(int(b), v, p))
            else:
                if not __silent: self.logger.info('  {1} ({2:.1f}%) contacts have {0} reads'.format(int(b), v, p))

        if outfile:
            oh = gzip.open(outfile, 'wt')
            oh.write('{0}\n'.format('\t'.join(['chrom1', 'left1', 'right1', 'chrom2', 'left1', 'right2', 'read_count'])))
            for contact in store:
                if store[contact] >= threshold:
                    oh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(contact[0], contact[1], contact[1]+window,
                        contact[2], contact[3], contact[3]+window,
                        store[contact]))
            oh.close()

            if not __silent: self.logger.info('Saved {0}'.format(outfile))

        # output for contact_zscore_cov
        return dict(zip([int(i) for i in h[1]], h[0])), peaklen, len_peaks

    def bed_to_genespromoter(self,
        reads,
        bed,
        genome_data,
        window,
        outfile,
        threshold=1, # in reads;
        **kargs):
        '''
        Measure the contacts between a BED file and promoter;

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
                peaks[chrom] = set([])

            # We bin the peak to the nearest window
            cpt = (int(peak[1]) + int(peak[2])) // 2
            loc = (cpt // window) * window

            peaks[chrom].add(loc)
        bedin.close()

        self.logger.info('Found {0:,} BED peaks'.format(len_peaks))

        # And do the same for the genes;
        genes = {}
        for len_genes, gene in enumerate(genome_data):
            chrom = 'chr{0}'.format(gene['loc'].loc['chr'])

            if chrom not in genes:
                genes[chrom] = set([])

            # We bin the peak to the nearest window
            cpt = gene['loc'].loc['left'] # always a point;
            loc = (cpt // window) * window

            genes[chrom].add(loc)

        chrom_intersect = set(genes.keys()) & set(peaks.keys())

        for chrom in chrom_intersect:
            # Fill in if chrom not present in other list;;
            if chrom not in genes:
                continue
            if chrom not in peaks:
                continue

            genes_chrom = genes[chrom]
            peaks_chrom = peaks[chrom]

            common = genes_chrom.intersection(peaks_chrom)
            peaks[chrom] = peaks_chrom - common
            genes[chrom] = genes_chrom - common

        self.logger.info('Shrunk peaks from {0:,} to {1:,} by removing duplicate bins and common bins'.format(len_peaks, sum([len(peaks[chrom]) for chrom in peaks])))
        self.logger.info('Shrunk genes from {0:,} to {1:,} by removing duplicate bins and common bins'.format(len_genes, sum([len(genes[chrom]) for chrom in genes])))

        store = {}

        for idx, read_pair in enumerate(readsin):
            if (idx+1) % 1e6 == 0:
                self.logger.info('{0:,} reads processed'.format(idx+1))
                #break

            read_pair = read_pair.strip().split('\t')

            chrom_left = read_pair[0]
            chrom_rite = read_pair[3]

            cpt = (int(read_pair[1]) + int(read_pair[2])) // 2
            bin_left = (cpt // window) * window
            cpt = (int(read_pair[4]) + int(read_pair[5])) // 2
            bin_rite = (cpt // window) * window

            try:
                # See if there is a contact between peaks and genes:
                if bin_left in peaks[chrom_left] and bin_rite in genes[chrom_rite]:
                    korder = (chrom_left, bin_left, chrom_rite, bin_rite)

                elif bin_rite in peaks[chrom_rite] and bin_left in genes[chrom_left]:
                    korder = (chrom_rite, bin_rite, chrom_left, bin_left)

                else: # no contact;
                    continue

                if korder not in store:
                    store[korder] = 0
                store[korder] += 1

            except KeyError:
                pass # chrom is not in peaks, but is in reads

        readsin.close()

        # Work out histogram;
        hist_max = 21
        all_scores = list(store.values())
        h = numpy.histogram(all_scores, range=[1,hist_max], bins=hist_max-1)
        self.logger.info('Histogram of contacts:')
        tot = sum(i for i in h[0])
        perc = [i/tot*100 for i in h[0]]
        for v, b, p in zip(h[0], h[1], perc):
            if int(b) == hist_max-1:
                if int(b) == 1:
                    self.logger.info('  {1} ({2:.1f}%) contacts have {0}+ read'.format(int(b), v, p))
                else:
                    self.logger.info('  {1} ({2:.1f}%) contacts have {0}+ reads'.format(int(b), v, p))
            else:
                self.logger.info('  {1} ({2:.1f}%) contacts have {0} reads'.format(int(b), v, p))

        oh = gzip.open(outfile, 'wt')
        oh.write('{0}\n'.format('\t'.join(['bedchrom1', 'bedleft1', 'bedright1', 'promchrom2', 'promleft1', 'promright2', 'read_count'])))
        for contact in store:
            if store[contact] >= threshold:
                oh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(contact[0], contact[1], contact[1]+window,
                    contact[2], contact[3], contact[3]+window,
                    store[contact]))

        self.logger.info('Saved {0}'.format(outfile))

        oh.close()

        return
