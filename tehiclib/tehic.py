
import os
import gzip

from . import common
from .collect_valid_pairs import collect_valid_pairs, save_valid_pairs
from .assign_to_te import map_pairs
from .quantify_links import quantify
from .build_matrices import build_matrices
from .genome import valid_assemblies

from . import miniglbase3

class te_hic:
    def __init__(self,
        genome:str,
        label:str,
        logger,
        save_intermediate_files:bool=False
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
        assert genome in valid_assemblies, f'{genome} not in one of the valid genomes {common.valid_assemblies}'

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

        self.valid_chroms = sorted(self.chrom_sizes.keys())
        self.logger.info(f'valid chromosome are: {", ".join(self.valid_chroms)}')
        self.valid_chroms = set(self.valid_chroms)

        # Check they exist, but don't load until stage2 to save memory
        self.genome = os.path.join(self.__script_path, f'../genome/{genome}_tes.glb')
        assert os.path.exists(self.genome), f'{genome} not found'

        # Get the TE frequencies tables
        self.te_genome_freqs = miniglbase3.glload(os.path.join(self.__script_path, f'../genome/{genome}_te_genome_freqs.glb'))

        return

    def stage1_collect_valid_pairs(self, bam1_filename, bam2_filename, min_dist=5000, min_qual=10):
        '''
        **Stage 1**
            Collect QC passing valid pairs
        '''

        self.valid_pairs_tmp_file = collect_valid_pairs(
            bam1_filename, bam2_filename,
            min_dist=5000, min_qual=min_qual,
            label=self.label,
            logger=self.logger,
            valid_chroms=self.valid_chroms,
            _save_intermediate_files=self.__save_intermediate_files
            )

        return True

    def __set_to_str(self, set_obj):
        # It's already a str at this point, not a real set
        return str(set_obj).translate({ord(c): None for c in "{}, '"}).replace('set()', 'None')

    def stage2_assign_to_genome_feature(self):
        '''
        **Stage 2**
            Assign reads to a genome feature
        '''
        assert self.valid_pairs_tmp_file, 'Stage 1 results "valid pairs" has not been generated correctly'

        self.genome = miniglbase3.glload(os.path.join(self.__script_path, self.genome))

        self.mapped_pairs_temp_file = map_pairs(self.valid_pairs_tmp_file, genome=self.genome, label=self.label, logger=self.logger)

        del self.genome

        # TODO: convert the temp to a BEDPE;
        # save out the BEDPEs:
        self.logger.info('Saving BEDPE files')
        file_all = gzip.open(f'stage2.all.{self.label}.bedpe.gz', 'wt')
        file_tete = gzip.open(f'stage2.tete.{self.label}.bedpe.gz', 'wt')
        file_tenn = gzip.open(f'stage2.tenn.{self.label}.bedpe.gz', 'wt')
        file_nnnn = gzip.open(f'stage2.nnnn.{self.label}.bedpe.gz', 'wt')
        mapped_pairs_temp_file = open(self.mapped_pairs_temp_file, 'r')

        # File format:
        # (chromA, leftA, riteA, chromB, leftB, riteB, read1_feat, read1_type, read2_feat, read2_type, read1_strand, read2_strand)

        for done, pair in enumerate(mapped_pairs_temp_file):
            pair = pair.strip().split('\t')

            line = f'chr{pair[0]}\t{pair[1]}\t{pair[2]}\tchr{pair[3]}\t{pair[4]}\t{pair[5]}\t0\t{pair[10]}\t{pair[11]}\t{self.__set_to_str(pair[7])}-{self.__set_to_str(pair[6])}\t{self.__set_to_str(pair[9])}-{self.__set_to_str(pair[8])}\n'

            # All:
            file_all.write(line)

            if 'TE' in pair[7] and 'TE' in pair[9]: # TE <=> TE
                file_tete.write(line)
            elif 'TE' in pair[7] or 'TE' in pair[9]: # TE <=> non-TE
                file_tenn.write(line)
            else: # non-TE <=> non-TE
                file_nnnn.write(line)

        mapped_pairs_temp_file.close()
        file_all.close()
        file_tete.close()
        file_tenn.close()
        file_nnnn.close()

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
