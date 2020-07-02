#!/usr/bin/env python3
'''

Pair up all valid reads from a pair of bams

'''

import sys, os, gzip
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

    **Returns**
        A list of valid fastq pairs
    '''
    assert bam1_filename, 'You must provide a valid filename in bam1_filename'
    assert bam2_filename, 'You must provide a valid filename in bam2_filename'

    print('Started %s and %s' % (bam1_filename, bam2_filename))

    bf1 = pysam.AlignmentFile(bam1_filename, 'rb')
    bf2 = pysam.AlignmentFile(bam2_filename, 'rb')
    pairs = set([])
    pairs_add = pairs.add # speedup to skip binding

    # We assume the bam files are sorted by name and unaligned were also output
    stats_total_reads = 0
    stats_lowq = 0
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
            print('ERROR: Mismatched read names (%s != %s), make sure the BAMs contain unaligned reads and are sorted by name' % (read1.query_name, read2.query_name))
            sys.quit()

        # check both reads are aligned
        if read1.is_unmapped and read2.is_unmapped:
            stats_unaligned += 1
            continue
        elif read1.is_unmapped or read2.is_unmapped:
            stats_1aligned += 1
            continue
        stats_aligned += 1

        # Check both reads have reasonable quality
        if read1.mapping_quality < 10 and read2.mapping_quality < 10:
            stats_lowq += 1
            continue

        # Now trim if not a valid chrom (so the above stats are correct)
        if read1.reference_name not in valid_chroms:
            continue
        if read2.reference_name not in valid_chroms:
            continue

        # criteria1: Check the distance between the two reads
        dist = max([abs(read1.reference_start - read2.reference_end), abs(read1.reference_end - read2.reference_start)])
        if dist < min_dist:
            reject_too_close += 1
            continue

        # This does duplicate removal in one go.
        # Only use the starts, it saves memory, and anyway the 3' ends are unreliable for duplicate removal if the input has been quality/adapter trimmed
        pairs_add((read1.reference_name, read1.reference_start, read2.reference_name, read2.reference_end))
        done += 1 # subtract this number to get the number of duplicates removed
        #if done > 200000:
        #    break

        if stats_total_reads % 1000000 == 0:
            print('Processed: {:,}'.format(stats_total_reads))

    bf1.close()
    bf2.close()

    print('\ncollect_valid_pairs() stats:')
    print('  Aligned:')
    print('    Reads processed : {:,}'.format(stats_total_reads))
    print('    Low quality     : {:,} ({:.2%})'.format(stats_lowq, stats_lowq/stats_total_reads))
    print('    Correctly paired: {:,} ({:.2%})'.format(stats_aligned, stats_aligned/stats_total_reads))
    print('    One pair aligned: {:,} ({:.2%})'.format(stats_1aligned, stats_1aligned/stats_total_reads))
    print('    No pairs aligned: {:,} ({:.2%})'.format(stats_unaligned, stats_unaligned/stats_total_reads))
    print('  Criteria rejected:')
    print('    Too close       : {:,} ({:.2%})'.format(reject_too_close, reject_too_close/stats_total_reads))
    print('    Duplicates      : {:,} ({:.2%})'.format(done-len(pairs), (done-len(pairs))/done))
    print('  Final:')
    print('    Kept reads      : {:,} ({:.2%})'.format(len(pairs), len(pairs)/stats_total_reads))
    return pairs

def remove_duplicates(pairs):
    '''
    **Purpose**
        Remove exact duplicates, as likely PCR errors

        The naive algoritm (set(pairs))
        fails with out of memory if you are ~180M items...
    '''
    newp = set()
    for p in pairs:
        if p not in newp:
            newp.add(p)

    print('\nremove_duplicates() stats:')
    print('  Duplicates        : {:,} ({:.2%})'.format(len(pairs)-len(newp), (len(pairs)-len(newp))/len(pairs)))
    return list(newp)

def save_valid_pairs(pairs, output):
    '''
    **Purpose**
        Save the valid pairs to output
    '''
    oh = gzip.open(output, 'w')
    #oh.write('%s\n' % '\t'.join(['chrom1', 'start', 'end', 'chr2', 'start', 'end']))
    for p in pairs:
        oh.write('%s\n' % '\t'.join([p[0], str(p[1]), str(p[1]+50), p[2], str(p[3]-50), str(p[3])]))
    oh.close()
    print('    Saved           : {:,} pairs'.format(len(pairs)))

if __name__ == '__main__':
    work_path = sys.argv[0]

    if len(sys.argv) != 4:
        print('\nNot enough arguments!')
        print('collect_valid_pairs.py bam1 bam2 output.bedpe.gz')
        print('Also note the bam files MUST be sorted by name, or this will result in a mess')
        print()
        sys.exit()

    bam1_filename = sys.argv[1]
    bam2_filename = sys.argv[2]
    output = sys.argv[3]

    pairs = collect_valid_pairs(bam1_filename, bam2_filename)
    #pairs = remove_duplicates(pairs)
    save_valid_pairs(pairs, output)


