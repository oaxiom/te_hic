
import gzip
import sys
import os

from glbase3 import genome

if len(sys.argv) != 1:
    print('measure_gc_in_bed.py BEDFILE')

hg38 = genome()
hg38.bindSequence(os.path.expanduser('~/hg38/seq/'), memorymap=True)

filename = sys.argv[1]

if filename.endswith('.gz'):
    output = gzip.open(filename.replace('.bed', '.gc.bed'), 'wt')
    input = gzip.open(filename, 'rt')
else:
    output = open(filename.replace('.bed', '.gc.bed'), 'wt')
    input = open(filename, 'rt')

print(f"Processing {filename} to {filename.replace('.bed', '.gc.bed')}")

for line in input:
    line = line.strip().split('\t')

    chrom = line[0]
    l = int(line[1])
    r = int(line[2])

    seq = hg38.getSequence(f'chr{chrom}:{l}-{r}').upper()

    perc_gc = (seq.count('G') + seq.count('C')) / (r -l)
    perc_gc = int(perc_gc * 10)
    perc_gc *= 10

    output.write(f'{chrom}\t{l}\t{r}\t-\t{perc_gc}\n')

input.close()
output.close()