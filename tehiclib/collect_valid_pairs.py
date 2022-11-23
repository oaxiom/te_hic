'''

Pair up all valid reads from a pair of bams

'''

import sys, os, gzip
import pysam

valid_chroms = set(['chrX', 'chrY'] + [f'chr{i}' for i in range(1, 30)]) # cut scaffolds

def collect_valid_pairs(bam1_filename, bam2_filename, min_dist=5000, logger=None, min_qual=None):
    '''
    **Purpose**
        Collect valid pairs for some set of criteria

    **Arguments**
        bam1_filename, bam1_filename (Required)
            bam filenames to open

        Criteria implemented:
            min_dist (default=5000)
                minimum distance for a valid pair

            min_qual (default=10)
                minimum quality score

    **Returns**
        A list of valid fastq pairs
    '''
    assert bam1_filename, 'You must provide a valid filename in bam1_filename'
    assert bam2_filename, 'You must provide a valid filename in bam2_filename'
    assert min_qual, 'You must specify a minimum quality score'

    logger.info(f'Started {bam1_filename} and {bam2_filename}')

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
    stats_short_range = 0
    stats_long_range = 0
    stats_not_canonical_chromosome = 0
    reject_diff_chrom = 0
    reject_too_close = 0

    done = 0
    for read1, read2 in zip(bf1, bf2): # needs to be eof...
        stats_total_reads += 1
        # read name sanity check:
        if read1.query_name != read2.query_name:
            logger.error('Mismatched read names (%s != %s), make sure the BAMs contain unaligned reads and are sorted by name' % (read1.query_name, read2.query_name))
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
        if read1.mapping_quality < min_qual and read2.mapping_quality < min_qual:
            stats_lowq += 1
            continue

        # Now trim if not a valid chrom (so the above stats are correct)
        if read1.reference_name not in valid_chroms:
            stats_not_canonical_chromosome += 1
            continue
        if read2.reference_name not in valid_chroms:
            stats_not_canonical_chromosome += 1
            continue

        # criteria1: Check the distance between the two reads
        dist = max([abs(read1.reference_start - read2.reference_end), abs(read1.reference_end - read2.reference_start)])
        if dist < min_dist:
            reject_too_close += 1
            continue

        if dist < 20000:
            stats_short_range += 1
        elif dist > 20000:
            stats_long_range += 1

        loc_strand1 = '-' if read1.is_reverse else '+'
        loc_strand2 = '-' if read2.is_reverse else '+'

        # This does duplicate removal in one go.
        # Only use the starts, it saves memory, and anyway the 3' ends are unreliable for duplicate removal if the input has been quality/adapter trimmed
        # Also strip the 'chr' off the front of the contigs. Could be a problem for some genomes?
        pairs_add((read1.reference_name[3:], read1.reference_start, read2.reference_name[3:], read2.reference_end, loc_strand1, loc_strand2))
        done += 1 # subtract this number to get the number of duplicates removed

        if stats_total_reads % 1000000 == 0:
            logger.info('Processed: {:,}'.format(stats_total_reads))

    bf1.close()
    bf2.close()

    logger.info('\ncollect_valid_pairs() stats:')
    logger.info('  Aligned:')
    logger.info('    Reads processed           : {:,}'.format(stats_total_reads))
    logger.info('    Correctly paired          : {:,} ({:.2%})'.format(stats_aligned, stats_aligned/stats_total_reads))
    logger.info('  Rejected reads:')
    logger.info('    Low quality               : {:,} ({:.2%})'.format(stats_lowq, stats_lowq/stats_total_reads))
    logger.info('    One pair aligned          : {:,} ({:.2%})'.format(stats_1aligned, stats_1aligned/stats_total_reads))
    logger.info('    Not canonical chromosome  : {:,} ({:.2%})'.format(stats_not_canonical_chromosome, stats_not_canonical_chromosome/stats_total_reads))
    logger.info('    No pairs aligned          : {:,} ({:.2%})'.format(stats_unaligned, stats_unaligned/stats_total_reads))
    logger.info('    Duplicates                : {:,} ({:.2%})'.format(done-len(pairs), (done-len(pairs))/done))
    logger.info('  Rejected reads (by criteria):')
    logger.info('    Too close                 : {:,} ({:.2%})'.format(reject_too_close, reject_too_close/stats_total_reads))
    logger.info('  Final:')
    logger.info('    Kept reads                : {:,} ({:.2%})'.format(len(pairs), len(pairs)/stats_total_reads))
    logger.info('    Kept short-range (<20kb)  : {:,} ({:.2%})'.format(stats_short_range, stats_short_range/stats_total_reads))
    logger.info('    Kept long-range (>20kb)   : {:,} ({:.2%})'.format(stats_long_range,  stats_long_range/stats_total_reads))

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
    print('  Duplicates                : {:,} ({:.2%})'.format(len(pairs)-len(newp), (len(pairs)-len(newp))/len(pairs)))
    return list(newp)

def save_valid_pairs(pairs, output):
    '''
    **Purpose**
        Save the valid pairs to output
    '''
    oh = gzip.open(output, 'wt')
    #oh.write('%s\n' % '\t'.join(['chrom1', 'start', 'end', 'chr2', 'start', 'end']))
    for p in pairs:
        oh.write('%s\n' % '\t'.join([p[0], str(p[1]), str(p[1]+50), p[2], str(p[3]-50), str(p[3]), '.', '0', p[4], p[5]]))
    oh.close()
