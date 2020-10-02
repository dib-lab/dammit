# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

import click
import pandas as pd

from ope.io.hmmer import HMMerParser
from ope.io.gff3 import GFF3Parser 
from ope.io.transdecoder import TransDecoderPepParser


def split_transdecoder_names(hmmer_df):
    hmmer_df.rename(columns={'query_name': 'full_query_name'},
                    inplace=True)

    name_df = hmmer_df.full_query_name.str.split('::', expand=True)
    name_df.rename(columns={0: 'td_transcript_name',
                            1: 'query_name',
                            2: 'td_g',
                            3: 'td_m'},
                   inplace=True)
    hmmer_df = pd.concat((hmmer_df, name_df), axis=1)
    
    return hmmer_df


@click.command('remap-hmmer-coords')
@click.argument('hmmerfilename')
@click.argument('pep_fa_filename')
@click.argument('output_filename')
def remap_hmmer_coords_cmd(hmmer_filename, pep_fa_filename, output_filename):
    '''Remap hmmscan coordinates using TransDecoder ORF coordinates.

    HMMER_FILENAME: HMMER results given by alignment to TransDecoder ORFS.
    PEP_FA_FILENAME: TransDecoder ORFS given by TransDecoder.LongOrfs.
    '''


    ref_df = TransDecoderPepParser(pep_fa_filename).read()
    hmmer_df = HMMerParser(hmmer_filename).read()

    if len(ref_df) > 0 and len(hmmer_df) > 0:
        merged_df = pd.merge(hmmer_df,
                                ref_df,
                                left_on='query_name',
                                right_on='full_transcript_name')

        hmmer_df['env_coord_from'] = (merged_df.src_start + \
                                        (3 * merged_df.env_coord_from)).astype(int)
        hmmer_df['env_coord_to'] = (merged_df.src_start + \
                                    (3 * merged_df.env_coord_to)).astype(int)
        hmmer_df['ali_coord_from'] = (merged_df.src_start + \
                                        (3 * merged_df.ali_coord_from)).astype(int)
        hmmer_df['ali_coord_to'] = (merged_df.src_start + \
                                    (3 * merged_df.ali_coord_to)).astype(int)
        hmmer_df['strand'] = merged_df['strand']
        hmmer_df['full_feature_name'] = merged_df['full_transcript_name']
        hmmer_df['query_name'] = merged_df['transcript_name']

    hmmer_df.to_csv(output_filename, header=True, index=False)
