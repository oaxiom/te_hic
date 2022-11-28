
import os, gzip
from . import common
from .collect_valid_pairs import collect_valid_pairs, save_valid_pairs
from .assign_to_te import map_pairs
from .quantify_links import quantify

from . import miniglbase3

class te_hic:
    def __init__(self,
        genome:str,
        label:str,
        logger,
        save_intermediate_files=True
        ):
        '''
        Entry point for te_hic

        **Arguments**
            genome (Required)
                One of the supported genomes

            label (Required)
                sample label for saving files

            logger (Required)
                python logger for console output

            save_intermediate_files (Optional, default=True)
                Save the intermediate files from each stage.
                Mainly for debugging purposes.

        '''
        assert genome in common.valid_assemblies, f'{genome} not in one of the valid genomes {common.valid_assemblies}'

        self.label = label
        self.valid_pairs = None
        self.logger = logger
        self.__save_intermediate_files = save_intermediate_files
        self.__script_path = os.path.dirname(os.path.realpath(__file__))

        self.bind_genome(os.path.join(self.__script_path, f'../genome/{genome}_glb_gencode_tes.glb'))

        return

    def bind_genome(self, genelist_glb_filename):
        self.genome = miniglbase3.glload(genelist_glb_filename)

    def stage1_collect_valid_pairs(self, bam1_filename, bam2_filename, min_dist=5000, min_qual=10):
        '''
        **Stage 1**
            Collect QC passing valid pairs
        '''

        self.valid_pairs = collect_valid_pairs(bam1_filename, bam2_filename,
            min_dist=5000, min_qual=min_qual,
            logger=self.logger)

        if self.__save_intermediate_files:
            save_valid_pairs(self.valid_pairs, f'stage1.int.{self.label}.tsv.gz')
            self.logger.info(f'Intermediate file: Saved {len(self.valid_pairs):,} pairs')

        return True

    def stage2_assign_to_genome_feature(self):
        '''
        **Stage 2**
            Assign reads to a genome feature
        '''
        assert self.genome, 'genome must be valid'
        assert self.valid_pairs, 'Stage 1 results "valid pairs" has not been generated correctly'

        self.mapped_pairs = map_pairs(self.valid_pairs, genome=self.genome)

        del self.valid_pairs # not needed

        if self.__save_intermediate_files:
            # out form:
            # ((pairs, read1_feat, read1_type, read2_feat, read2_type))
            # pairs = ('chr7', 150285954, 'chr4', 130529111, '-', '+');
            #
            out = gzip.open(f'stage2.int.{self.label}.tsv.gz', 'wt')
            for o in self.mapped_pairs:
                read1_names = ', '.join(o[1]) if o[1] else 'None'
                read1_types = ', '.join(o[2]) if o[2] else 'None'
                read2_names = ', '.join(o[3]) if o[3] else 'None'
                read2_types = ', '.join(o[4]) if o[4] else 'None'

                line = [f'chr{o[0][0]}', str(o[0][1]), str(o[0][1]+50),
                    read1_names, read1_types,
                    f'chr{o[0][2]}', str(o[0][3]-50), str(o[0][4]),
                    read2_names, read2_types]

                out.write('{}\n'.format('\t'.join(line)))
            out.close()
            self.logger.info(f'Intermediate file: Pair assignments {len(self.mapped_pairs):,} saved')

        return True

    def stage3_quantify_links(self):
        '''
        **Stage 3**
            Aggregate the genome contacts and output some statistics

        '''
        assert self.mapped_pairs, 'Stage 2 results "mapped_pairs" was not generated correctly'

        qfy = quantify(self.label, logger=self.logger)
        qfy.bind_genome(self.genome)
        qfy.measure_te_anchors(self.mapped_pairs)

        #del self.mapped_pairs
        return True

    def stage4_build_matrices(self):
        '''
        **Stage 4**
            Build the matrices at the required resolutions

        '''
        assert self.mapped_pairs, 'Stage 2 results "mapped_pairs" was not generated correctly'

        return True
