

import os
import glob
import pickle



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

def load_randoms(filename):
    for


if __name__ == '__main__':
    all_data = dict(
        reals = load_intercons('./real_con/*.intracon_num.txt'),
        bkgds = load_intercons('./random_con/*.intracon_num.txt'), # Liyang's background
        peaklens = load_peaklens('./cov_ctcf.txt'),
        randoms = load_randoms('./random/*.bed')
    )

    with open('./all_data.pkl', 'wb') as oh:
        pickle.dump(all_data, oh)

    print(all_data)

