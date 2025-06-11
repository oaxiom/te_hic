

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

def load_randoms(files):
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

if __name__ == '__main__':
    all_data = dict(
        reals = load_intercons('./real_con/*.intracon_num.txt'),
        bkgds = load_intercons('./random_con/*.intracon_num.txt'), # Liyang's background
        peaklens = load_peaklens('./cov_ctcf.txt'),
        randoms = load_randoms('./randoms/*.bed.gz')
    )

    with open('./all_data.pkl', 'wb') as oh:
        pickle.dump(all_data, oh)

    #print(all_data)

