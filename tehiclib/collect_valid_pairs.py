'''

Pair up all valid reads from a pair of bams

'''

import sys, os, gzip, random, subprocess, time
import pysam

valid_chroms = set(['chrX', 'chrY'] + [f'chr{i}' for i in range(1, 30)]) # cut scaffolds

def dump_to_file(pairs, filehandle):
    [filehandle.write(p) for p in sorted(pairs)]
    num_saved = len(pairs)
    del pairs
    return num_saved

def collect_valid_pairs(bam1_filename,
    bam2_filename,
    min_dist=5000,
    label=None,
    logger=None,
    _save_intermediate_files=False,
    min_qual=None):
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
    assert label, 'Need a label'

    logger.info(f'Started {bam1_filename} and {bam2_filename}')

    bf1 = pysam.AlignmentFile(bam1_filename, 'rb')
    bf2 = pysam.AlignmentFile(bam2_filename, 'rb')

    temp_filename = f'stage1.{str(int(time.time()))[5:]}{random.randint(10, 1e6):0>7}.{label}.tmp'
    temp_out = open(temp_filename, 'w')

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

    pairs = set([])
    pairs_add = pairs.add

    step = int(100e6) # Peak memory ~10Gb
    total_saved = 0

    done = 0
    for read1, read2 in zip(bf1, bf2): # needs to be eof...
        stats_total_reads += 1
        if stats_total_reads % step == 0:
            num_saved = dump_to_file(pairs, temp_out) # semi pair removed;
            pairs = set([])
            pairs_add = pairs.add
            total_saved += num_saved # reset so last write is accurate
            logger.info(f'Processed: {stats_total_reads:,} reads, removed {step-num_saved:,} ({(step-num_saved)/step:.1%}) reads by preduplicate removal')

        # read name sanity check:
        if read1.query_name != read2.query_name:
            logger.error(f'Mismatched read names ({read1.query_name} != {read2.query_name}), make sure the BAMs contain unaligned reads and are sorted by name')
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
        #temp_out.write(f'{read1.reference_name[3:]}\t{read1.reference_start}\t{read2.reference_name[3:]}\t{read2.reference_end}\t{loc_strand1}\t{loc_strand2}\n')
        pairs_add(f'{read1.reference_name[3:]}\t{read1.reference_start}\t{read2.reference_name[3:]}\t{read2.reference_end}\t{loc_strand1}\t{loc_strand2}\n')
        done += 1

    bf1.close()
    bf2.close()

    # Final dump
    num_saved = dump_to_file(pairs, temp_out)
    total_saved += num_saved
    logger.info(f'Processed: {stats_total_reads:,} reads, removed {stats_total_reads-total_saved:,} ({(stats_total_reads-total_saved)/stats_total_reads:.1%}) reads by preduplicate removal')
    temp_out.close()
    del pairs

    logger.info('Sorting temp file')

    # It seems awk gets killed in the region of ~200M so just use sort | uniq
    subprocess.run(f'sort {temp_filename} | uniq > {temp_filename}.sorted', shell=True)    # Portable memory resilient, slower?
    #if done > 200e6: # We saw failures at 650M and 1300M
    #subprocess.run(f"awk '!x[$0]++' {temp_filename} > {temp_filename}.sorted", shell=True) # Faster, higher peak memory?!

    if not _save_intermediate_files:
        os.remove(f'{temp_filename}') # only sorted needed now;

    #logger.info('Reloading')
    #pairs = set([])
    #pairs_add = pairs.add # speedup to skip binding
    #oh = open(f'{temp_filename}.sorted', 'r')
    #for line in oh:
    #    line = line.strip().split('\t')
    #    pairs_add((line[0], int(line[1]), line[2], int(line[3]))) # You can drop the strands now as not needed anymore, line[4], line[5]))
    #oh.close()

    logger.info('Stage 1 stats:')
    logger.info('  Aligned:')
    logger.info('    Reads processed           : {:,}'.format(stats_total_reads))
    logger.info('    Correctly paired          : {:,} ({:.2%})'.format(stats_aligned, stats_aligned/stats_total_reads))
    logger.info('  Rejected reads:')
    logger.info('    Low quality               : {:,} ({:.2%})'.format(stats_lowq, stats_lowq/stats_total_reads))
    logger.info('    One pair aligned          : {:,} ({:.2%})'.format(stats_1aligned, stats_1aligned/stats_total_reads))
    logger.info('    Not canonical chromosome  : {:,} ({:.2%})'.format(stats_not_canonical_chromosome, stats_not_canonical_chromosome/stats_total_reads))
    logger.info('    No pairs aligned          : {:,} ({:.2%})'.format(stats_unaligned, stats_unaligned/stats_total_reads))
    logger.info('  Rejected reads (by criteria):')
    logger.info('    Too close                 : {:,} ({:.2%})'.format(reject_too_close, reject_too_close/stats_total_reads))
    logger.info('  Final:')
    logger.info('    Kept reads                : {:,} ({:.2%})'.format(stats_total_reads-total_saved, (stats_total_reads-total_saved)/stats_total_reads))
    logger.info('    [Note that these numbers below include PCR duplicates]')
    logger.info('    Kept short-range (<20kb)  : {:,} ({:.2%})'.format(stats_short_range, stats_short_range/stats_total_reads))
    logger.info('    Kept long-range (>20kb)   : {:,} ({:.2%})'.format(stats_long_range,  stats_long_range/stats_total_reads))

    ret = subprocess.run(f"wc {temp_filename}.sorted", shell=True, capture_output=True).stdout.decode()
    ret = int(str(ret).strip().split(' ')[0])
    logger.info('    [Sorted, unique should match above number]')
    logger.info('    Kept long-range (>20kb)   : {:,} ({:.2%})'.format(ret,  ret/stats_total_reads))

    return f"{temp_filename}.sorted"

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
