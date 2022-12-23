
import os, gzip
from . import common
from .collect_valid_pairs import collect_valid_pairs, save_valid_pairs
from .assign_to_te import map_pairs
from .quantify_links import quantify
from .build_matrices import build_matrices

from . import miniglbase3

class te_hic:
    def __init__(self,
        genome:str,
        label:str,
        logger,
        save_intermediate_files=False
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
        self.valid_pairs_tmp_file = None
        self.logger = logger
        self.__save_intermediate_files = save_intermediate_files
        self.__script_path = os.path.dirname(os.path.realpath(__file__))


        oh = open(os.path.join(self.__script_path, f'../genome/{genome}.chromSizes.clean'), 'r')
        self.chrom_sizes = {}
        for line in oh:
            line = line.strip().split('\t')
            self.chrom_sizes[line[0]] = int(line[1])
        oh.close()

        # Check they exist, but don't load until stage2 to save memory
        self.genome = os.path.join(self.__script_path,f'../genome/{genome}_glb_gencode_tes.glb')
        assert os.path.exists(self.genome), f'{genome} not found'

        # Get the TE frequencies tables
        self.te_genome_freqs = miniglbase3.glload(os.path.join(self.__script_path, f'../genome/{genome}_te_genome_freqs.glb'))

        return

    def stage1_collect_valid_pairs(self, bam1_filename, bam2_filename, min_dist=5000, min_qual=10):
        '''
        **Stage 1**
            Collect QC passing valid pairs
        '''

        self.valid_pairs_tmp_file = collect_valid_pairs(bam1_filename, bam2_filename,
            min_dist=5000, min_qual=min_qual,
            label=self.label,
            logger=self.logger,
            _save_intermediate_files=self.__save_intermediate_files)

        return True

    def stage2_assign_to_genome_feature(self):
        '''
        **Stage 2**
            Assign reads to a genome feature
        '''
        self.genome = miniglbase3.glload(os.path.join(self.__script_path, self.genome))

        assert self.valid_pairs_tmp_file, 'Stage 1 results "valid pairs" has not been generated correctly'

        self.mapped_pairs_temp_file = map_pairs(self.valid_pairs_tmp_file, genome=self.genome, label=self.label, logger=self.logger)

        if not self.__save_intermediate_files:
            os.remove(self.valid_pairs_tmp_file) # not needed anymore

        return True

    def stage3_quantify_links(self):
        '''
        **Stage 3**
            Aggregate the genome contacts and output some statistics

        '''
        assert self.mapped_pairs_temp_file, 'Stage 2 results "mapped_pairs" was not generated correctly'

        qfy = quantify(self.label, logger=self.logger)
        qfy.bind_te_freqs(self.te_genome_freqs)
        qfy.measure_te_anchors(self.mapped_pairs_temp_file)

        return True

    def stage4_build_matrices(self, resolutions):
        '''
        **Stage 4**
            Build the matrices at the required resolutions

        '''
        assert self.mapped_pairs_temp_file, 'Stage 2 results "mapped_pairs" was not generated correctly'

        for resolution in resolutions:
            mat = build_matrices(self.chrom_sizes, resolution, logger=self.logger)

            mat.build_matrices(self.mapped_pairs_temp_file)
            mat.save_matrices(self.label)

            del mat

        if not self.__save_intermediate_files:
            os.remove(self.mapped_pairs_temp_file) # not needed anymore

        return True
