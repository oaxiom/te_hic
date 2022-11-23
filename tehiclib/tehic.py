
import os
from . import common
from .collect_valid_pairs import collect_valid_pairs, remove_duplicates, save_valid_pairs

class te_hic:
    def __init__(self,
        label:str,
        logger,
        save_intermediate_files=True
        ):
        '''
        Entry point for te_hic

        **Arguments**
            label (Required)
                sample label for saving files

            logger (Required)
                python logger for console output

            save_intermediate_files (Optional, default=True)
                Save the intermediate files from each stage.
                Mainly for debugging purposes.

        '''
        self.label = label
        self.valid_pairs = None
        self.logger = logger
        self.__save_intermediate_files = save_intermediate_files
        self.__script_path = os.path.dirname(os.path.realpath(__file__))

        return

    def stage1_collect_valid_pairs(self, bam1_filename, bam2_filename, min_dist=5000, min_qual=10):
        '''
        **Stage 1**
            Collect QC passing valid pairs
        '''

        self.valid_pairs = collect_valid_pairs(bam1_filename, bam2_filename,
            min_dist=5000, min_qual=min_qual,
            logger=self.logger)

        if self.__save_intermediate_files:
            save_valid_pairs(self.valid_pairs, f'stage1.int.{self.label}.tsv')
            self.logger.info(f'Intermediate file: Saved {len(self.valid_pairs):,} pairs')

        return True

    def stage2_assign_to_genome_feature(self, species):
        '''
        **Stage 2**
            Assign reads to a genome feature
        '''
        assert species, 'species must be valid'
        assert species in common.valid_assemblies, f'species {species} not found in supported species {common.valid_assemblies}'
        if not self.valid_pairs:
            raise AssertionError('valid pairs has not been generated')

        mte = measureTE(species, logger=self.logger)
        mte.bind_genome(os.path.join(script_path, 'genome/%s_glb_gencode_tes.glb' % species))
        mte.load_bedpe(sys.argv[2], sys.argv[3])

        return True

    def stage3_quantify_links(self):
        '''
        **Stage 3**
            Aggregate the genome contacts.

        '''
        return True

    def stage4_build_matrices(self):
        '''
        **Stage 4**
            Build the matrices at the required resolutions

        '''
        return True
