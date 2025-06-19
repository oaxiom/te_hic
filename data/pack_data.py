

import os
import glob
import pickle
import gzip

def open_intracon_file(filename):
    cons = []
    with open(filename, 'r') as oh:
        for line in oh:
            if '#' in line:
                continue
            line = line.strip().split('\t')
            cons.append(int(line[1])) # 3 -> 21
    return cons

def load_intercons(path):
    con = {}
    for f in glob.glob(path):
        # I use 1 random here, whereas Liyang used 10 randoms in the paper.
        name = os.path.split(f)[1].replace('hesc_primed_', '').replace('.intracon_num.txt', '').split('.')[0]
        dat = open_intracon_file(f)
        con[name] = dat

    # Stuff them all into a big table;
    #contacts = pd.concat(con, axis=1)
    #contacts.columns = [ os.path.split(f)[1].replace('hesc_primed_','').replace('.intracon_num.txt','') for f in real_file ]

    return con

def load_peaklens(filename):
    ## load total length
    peaklens = {}
    with open(filename, 'r') as oh:
        for line in oh:
            if line.startswith('tf'):
                continue
            line = line.strip().split(' ')
            peaklens[line[0]] = int(line[2])

    #peaklen = dict(pd.read_csv('./cov_ctcf.txt', header=0, sep=' ')[['tf','len']].values)
    return peaklens

def load_beds(files):
    randoms = {}
    loci_loaded = 0
    for file in glob.glob(files):
        with gzip.open(file, 'rt') as oh:
            for line in oh:
                line = line.strip().split('\t')
                chrom = line[0]
                if chrom not in randoms:
                    randoms[chrom] = []
                randoms[chrom].append(int(line[1])) # I only need one point;
                loci_loaded += 1
    print(f'Found {loci_loaded} random loci')
    return randoms

def load_bed(file):
    randoms = {}
    loci_loaded = 0
    with gzip.open(file, 'rt') as oh:
        for line in oh:
            line = line.strip().split('\t')
            chrom = line[0]
            if chrom not in randoms:
                randoms[chrom] = []
            randoms[chrom].append(int(line[1])) # I only need one point;
            loci_loaded += 1

        for chrom in randoms: # Save space;
            randoms[chrom] = list(set(randoms[chrom]))

    print(f'Found {loci_loaded} random loci')
    return randoms

def load_randoms_gc(files, fasta):
    from glbase3 import genome # Need to manipulate FASTA;

    peak_size = 200

    hg38 = genome()
    hg38.bindSequence(fasta, memorymap=True)

    randoms_by_gc_percent = {20: {}, 30: {}, 40: {}, 50: {}, 60: {}, 70: {}}

    loci_loaded = 0
    for file in glob.glob(files):
        with gzip.open(file, 'rt') as oh:
            for line in oh:
                line = line.strip().split('\t')
                chrom = line[0]

                l = int(line[1])
                r = l + peak_size # ~Average peak size for narrow peaks;

                seq = hg38.getSequence(f'chr{chrom}:{l}-{r}').upper()
                perc_gc = (seq.count('G') + seq.count('C')) / peak_size
                perc_gc = int(perc_gc * 10)
                perc_gc *= 10

                # There are too few loci to get full coverage at these extreme GC contents, so bracket them.

                if perc_gc < 20: perc_gc = 20
                if perc_gc > 70: perc_gc = 70

                #if perc_gc == 10:
                #    print(perc_gc, chrom, l)
                #    print(seq)

                if chrom not in randoms_by_gc_percent[perc_gc]:
                    randoms_by_gc_percent[perc_gc][chrom] = []
                randoms_by_gc_percent[perc_gc][chrom].append(int(line[1])) # I only need one point;
                loci_loaded += 1
                if (loci_loaded+1) % 1e4 == 0:
                    print('{0:,} peaks processed'.format(loci_loaded+1))

    for perc in randoms_by_gc_percent:
        num_loci_this_perc = 0
        for chrom in randoms_by_gc_percent[perc]:
            num_loci_this_perc += len(randoms_by_gc_percent[perc][chrom])
        print(f'Found {num_loci_this_perc} random loci with GC% {perc}%')

    return randoms_by_gc_percent

if __name__ == '__main__':
    all_data = dict(
        peaklens=load_peaklens('./cov_ctcf.txt'),

        reals = load_intercons('./real_con/*.intracon_num.txt'),

        # This is how it was first done, with bedtools shuf
        bkgds = load_intercons('./random_con/*.intracon_num.txt'), # Liyang used 10 backgrounds per real, but 1 is enough.

        # For dynamic random generation
        randoms = load_beds('./randoms/*.bed.gz'),
        randoms_gc = load_randoms_gc('./randoms/*.bed.gz', os.path.expanduser('~/hg38/seq/')),
        randoms_pooled = load_bed('./peaks/all_peaks.bed.gz'), # A third normalisation technique, this time using pools of all peaks;

        # Intracons generated using matched GC bacgrounds
        bkgds_gc = load_intercons('./gc_random_con/*.intracon_num.txt'),
        #bkgds_pooled = load_intercons('./pooled_random_con/*.intracon_num.txt'),
    )

    with open('./all_data.pkl', 'wb') as oh:
        pickle.dump(all_data, oh)

    #print(all_data)

