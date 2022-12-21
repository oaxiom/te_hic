'''

Assign the BEDPE to a TE.

Does one or more end overlap with a TE?

'''

import sys, os, gzip, random
from . import miniglbase3
from . import common

def save_output(output_data, filename_handle):
    for line in output_data:
        filename_handle.write('{}\n'.format('\t'.join([str(i) for i in line])))

def map_pairs(valid_pairs_temp_file, genome, label=None, logger=False):
    '''
    **Purpose**
        Load in a BEDPE file, ideally output by collect_valid_pairs.py, although I guess any valid BEDPE will do

    **Arguments**
        valid_pairs (Required)
            A list containing valid pairs;

        genome (Required)
            genome glb file
    '''
    assert valid_pairs_temp_file, 'No valid pairs'
    assert logger, 'Need logger for output'
    assert label, 'Need a label'

    done = 0
    bucket_size = miniglbase3.config.bucket_size

    output = []

    self_genome_linearData = genome.linearData
    self_genome_buckets = genome.buckets

    #print(self_genome_buckets)

    valid_pairs = open(valid_pairs_temp_file, 'r')

    output_filename = f'stage2.{random.randint(10, 1e6):0>7}.{label}.tmp'
    output_file = open(output_filename, 'w')

    for idx, pairs in enumerate(valid_pairs):
        # pairs format ('chr7', 150285954, 'chr4', 130529111, '-', '+');
        pairs = pairs.strip().split('\t')
        chromA = pairs[0]
        leftA = int(pairs[1])
        riteA = leftA + 100

        left_buck = ((leftA-1)//bucket_size) * bucket_size
        right_buck = (riteA//bucket_size) * bucket_size
        buckets_reqd = range(left_buck, right_buck+bucket_size, bucket_size)
        result = []

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in self_genome_buckets[chromA]:
                    loc_ids.update(self_genome_buckets[chromA][buck]) # set = unique ids

            for index in loc_ids:
                if riteA >= self_genome_linearData[index]["loc"].loc["left"] and leftA <= self_genome_linearData[index]["loc"].loc["right"]:
                    result.append(self_genome_linearData[index])

            read1_feat = []
            read1_type = []
            if result:
                for r in result:
                    read1_feat.append(r['name'])
                    read1_type.append(r['type'])

        # work out which of the buckets is required:
        chromB = pairs[2]
        leftB = int(pairs[3]) - 100 # Yes, this is correct, it takes reference_end from the BAM
        riteB = int(pairs[3])

        left_buck = ((leftB-1)//bucket_size) * bucket_size
        right_buck = (riteB//bucket_size) * bucket_size
        buckets_reqd = range(left_buck, right_buck+bucket_size, bucket_size)
        result = []

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in self_genome_buckets[chromB]:
                    loc_ids.update(self_genome_buckets[chromB][buck]) # set = unique ids

            for index in loc_ids:
                if riteB >= self_genome_linearData[index]["loc"].loc["left"] and leftB <= self_genome_linearData[index]["loc"].loc["right"]:
                    result.append(self_genome_linearData[index])

            read2_feat = []
            read2_type = []
            if result:
                for r in result:
                    read2_feat.append(r['name'])
                    read2_type.append(r['type'])

        if read1_feat:
            read1_feat = set(read1_feat)
            read1_type = set(read1_type)
        else:
            read1_feat = set()
            read1_type = set()

        if read2_feat:
            read2_feat = set(read2_feat)
            read2_type = set(read2_type)
        else:
            read2_feat = set()
            read2_type = set()

        output.append((chromA, leftA, riteA, chromB, leftB, riteB, read1_feat, read1_type, read2_feat, read2_type))

        if idx % 1e6 == 0:
            logger.info(f'Processed: {idx:,}')
            # output reads results:
            save_output(output, output_file)
            output = []

    save_output(output, output_file)
    del output
    output_file.close()

    logger.info(f'Processed: {idx:,} reads in total')

    return output_filename
