#!/usr/bin/env python3
'''

Pair up all valid reads from a pair of bams

'''

import sys, os
import pysam

def collect_valid_pairs(bam1_filename, bam2_filename, min_dist=5000):
    '''
    **Purpose**
        Collect valid pairs for some set of criteria
        
    **Arguments**
        bam1_filename, bam1_filename (Required)
            bam filenames to open

        Criteria implemented:
            min_dist (default=5000)
                minimum distance for a valid pair

            inter_chrom (default=True)
                Only collect inter-chromosomes                
    
    **Returns**
        A list of valid fastq pairs
    '''
    assert bam1_filename, 'You must provide a valid filename in bam1_filename'
    assert bam2_filename, 'You must provide a valid filename in bam2_filename'

    bf1 = pysam.AlignmentFile(bam1_filename, 'rb')
    bf2 = pysam.AlignmentFile(bam2_filename, 'rb')
    output = []

    # We assume the bam files are sorted by name and unaligned were also output
    stats_aligned = 0
    stats_1aligned = 0
    stats_unaligned = 0
    reject_diff_chrom = 0 
    reject_too_close = 0

    done = 0
    while done < 10:
        read1 = bf1.__next__()
        read2 = bf2.__next__()    
        # read name sanity check:
        print(read1.query_name, read2.query_name)
        if read1.query_name != read2.query_name:
            print('ERORR: Mismatched read names (%s != %s), make sure the BAMs contain unaligned and are sorted by name' % (read1.query_name, read2.query_name))
            sys.quit()

        # First, check both reads are aligned
        if read1.is_unmapped and read2.is_unmapped:
            stats_unaligned += 1
            continue
        elif read1.is_unmapped or read2.is_unmapped:
            stats_1aligned += 1
            continue
        stas_aligned += 1
        
        # criteria1: Must be on same chrom
        if read1.reference_name != read2.reference_name:
            continue
    
        done += 1
        

    bf1.close()
    bf2.close()

    print('collect_valid_pairs() stats:')
    print('  Aligned:')
    print('    Reads processed : {:,}'.format(done))
    print('    Correctly paired: {:,}'.format(stats_aligned))
    print('    One pair aligned: {:,}'.format(stats_1aligned))
    print('    No pairs aligned: {:,}'.format(stats_unaligned))
    print('  Criteria rejected:')
    print('    Different chrom : {:,}'.format(reject_diff_chrom))
    print('    Too close       : {:,}'.format(reject_too_close))

    return output
    
def save_valid_pairs(pairs, output):
    pass

if __name__ == '__main__':
    work_path = sys.argv[0]

    if len(sys.argv) != 4:
        print('Not enough arguments!')
        print('2.pair.py bam1 bam2 output.bam')
        print('Also note the bam files MUST be sorted by name, or this will result in a mess')
        sys.exit()

    bam1_filename = sys.argv[1]
    bam2_filename = sys.argv[2]
    output = sys.argv[3]

    pairs = collect_valid_pairs(bam1_filename, bam2_filename)
    save_valid_pairs(pairs, output)


