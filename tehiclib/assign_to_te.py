#!/usr/bin/env python3

'''

Assign the BEDPE to a TE.

Does one or more end overlap with a TE?

'''

import sys, os, gzip
#import miniglbase3 # miniglbase3 namespace mangling!
from . import common

def map_pairs(valid_pairs, genome):
    '''
    **Purpose**
        Load in a BEDPE file, ideally output by collect_valid_pairs.py, although I guess any valid BEDPE will do

    **Arguments**
        valid_pairs (Required)
            A list containing valid pairs;

        genome (Required)
            genome glb file
    '''
    assert valid_pairs, 'No valid pairs'

    done = 0
    bucket_size = glbase3.config.bucket_size

    output = []

    self_genome_linearData = genome.linearData
    self_genome_buckets = genome.buckets

    for idx, pairs in enumerate(valid_pairs):
        # pairs format ('chr7', 150285954, 'chr4', 130529111, '-', '+');
        chrom = pairs[0]
        left = pairs[1]
        rite = rite + 50

        left_buck = ((left-1)//bucket_size) * bucket_size
        right_buck = (rite//bucket_size) * bucket_size
        buckets_reqd = range(left_buck, right_buck+bucket_size, bucket_size)
        result = []

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in self_genome_buckets[chrom]:
                    loc_ids.update(self_genome_buckets[chrom][buck]) # set = unique ids

            for index in loc_ids:
                if rite >= self_genome_linearData[index]["loc"].loc["left"] and left <= self_genome_linearData[index]["loc"].loc["right"]:
                    result.append(self_genome_linearData[index])

            read1_feat = []
            read1_type = []
            if result:
                for r in result:
                    read1_feat.append(r['name'])
                    read1_type.append(r['type'])

        # work out which of the buckets is required:
        chrom = pairs[2]
        left = pairs[3]
        rite = rite + 50

        left_buck = ((left-1)//bucket_size) * bucket_size
        right_buck = (rite//bucket_size) * bucket_size
        buckets_reqd = range(left_buck, right_buck+bucket_size, bucket_size)
        result = []

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in self_genome_buckets[chrom]:
                    loc_ids.update(self_genome_buckets[chrom][buck]) # set = unique ids

            for index in loc_ids:
                if rite >= self_genome_linearData[index]["loc"].loc["left"] and left <= self_genome_linearData[index]["loc"].loc["right"]:
                    result.append(self_genome_linearData[index])

            read2_feat = []
            read2_type = []
            if result:
                for r in result:
                    read2_feat.append(r['name'])
                    read2_type.append(r['type'])

        if read1_feat:
            read1_feat = ', '.join(set(read1_feat))
            read1_type = ', '.join(set(read1_type))
        else:
            read1_feat = 'None'
            read1_type = 'None'

        if read2_feat:
            read2_feat = ', '.join(set(read2_feat))
            read2_type = ', '.join(set(read2_type))
        else:
            read2_feat = 'None'
            read2_type = 'None'

        output.append((pairs, read1_feat, read1_type, read2_feat, read2_type))

        done += 1

        if done % 1000000 == 0:
            print(f'Processed: {done:,}')

    print(f'Processed {len(output):,} reads')

    return output
