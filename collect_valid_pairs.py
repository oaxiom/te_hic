#!/usr/bin/env python3
'''

Pair up all valid reads from a pair of bams

'''

import sys, os
import pysam

valid_chroms = set(['chrX', 'chrY'] + ['chr%s' % i for i in range(1, 30)]) # cut scaffolds

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

    print('Started %s and %s' % (bam1_filename, bam2_filename))

    bf1 = pysam.AlignmentFile(bam1_filename, 'rb')
    bf2 = pysam.AlignmentFile(bam2_filename, 'rb')
    pairs = []

    # We assume the bam files are sorted by name and unaligned were also output
    stats_total_reads = 0
    stats_aligned = 0
    stats_1aligned = 0
    stats_unaligned = 0
    stats_output = 0
    reject_diff_chrom = 0 
    reject_too_close = 0

    done = 0
    for read1, read2 in zip(bf1, bf2): # needs to be eof...
        stats_total_reads += 1
        # read name sanity check:
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
        stats_aligned += 1

        # Now trim if not a valid chrom (so the above stats are correct)
        if read1.reference_name not in valid_chroms:
            continue
        if read2.reference_name not in valid_chroms:
            continue
 
        # criteria1: Must be on same chrom
        #if read1.reference_name != read2.reference_name:
        #    reject_diff_chrom += 1
        #    continue
        # criteria2: Check the distance between the two reads
        dist = max([abs(read1.reference_start - read2.reference_end), abs(read1.reference_end - read2.reference_start)])
        if dist < min_dist:
            reject_too_close += 1
            continue

        pairs.append((read1.reference_name, read1.reference_start, read1.reference_end, read2.reference_name, read2.reference_start, read2.reference_end))  
    
        done += 1
        #if done > 200000:
        #    break

        if done % 1000000 == 0:
            print('Processed: {:,}'.format(done))
            
    bf1.close()
    bf2.close()

    print('\ncollect_valid_pairs() stats:')
    print('  Aligned:')
    print('    Reads processed : {:,}'.format(stats_total_reads))
    print('    Correctly paired: {:,} ({:.2%})'.format(stats_aligned, stats_aligned/stats_total_reads))
    print('    One pair aligned: {:,} ({:.2%})'.format(stats_1aligned, stats_1aligned/stats_total_reads))
    print('    No pairs aligned: {:,} ({:.2%})'.format(stats_unaligned, stats_unaligned/stats_total_reads))
    print('  Criteria rejected:')
    #print('    Different chrom : {:,} ({:.2%})'.format(reject_diff_chrom, reject_diff_chrom/stats_total_reads))
    print('    Too close       : {:,} ({:.2%})'.format(reject_too_close, reject_too_close/stats_total_reads))
    print('  Final:')
    print('    Kept reads      : {:,} ({:.2%})'.format(len(pairs), len(pairs)/stats_total_reads))
    return pairs
   
def remove_duplicates(pairs):
    '''
    **Purpose**
        Remove exact duplicates, as likely PCR errors
    '''
    newp = set(pairs)
    print('\nremove_duplicates() stats:')
    print('  Duplicates        : {:,} ({:.2%})'.format(len(pairs)-len(newp), (len(pairs)-len(newp))/len(pairs)))
    return list(newp)

def save_valid_pairs(pairs, output):
    '''
    **Purpose** 
        Save the valid pairs to output   
    '''
    oh = open(output, 'w')
    #oh.write('%s\n' % '\t'.join(['chrom1', 'start', 'end', 'chr2', 'start', 'end']))
    for p in pairs:
        oh.write('%s\n' % '\t'.join([p[0], str(p[1]), str(p[2]), p[3], str(p[4]), str(p[5])]))
    oh.close()
    print('\nsave_valid_pairs():')
    print('  Finally, saved : {:,} pairs'.format(len(pairs)))

if __name__ == '__main__':
    work_path = sys.argv[0]

    if len(sys.argv) != 4:
        print('\nNot enough arguments!')
        print('2.pair.py bam1 bam2 output.bedpe')
        print('Also note the bam files MUST be sorted by name, or this will result in a mess')
        print()
        sys.exit()

    bam1_filename = sys.argv[1]
    bam2_filename = sys.argv[2]
    output = sys.argv[3]

    pairs = collect_valid_pairs(bam1_filename, bam2_filename)
    pairs = remove_duplicates(pairs)
    save_valid_pairs(pairs, output)


