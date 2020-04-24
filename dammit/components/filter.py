# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import click
from shmlast import hits

from ..fileio.maf import MafParser


@click.command('best-hits')
@click.argument('maf_fn')
@click.argument('output_fn')
def maf_best_hits_cmd(maf_fn, output_fn):
    ''' Filter query best-hits from a MAF file.
    
    Get the best hits from a lastal MAF file: that is,
    for each query sequence, choose its top-scoring alignment.
    OUTPUT_FN is in csv format.
    \f

    Args:
        maf_fn (str): Path to the MAF file.
        output_fn (str): Path to store resulting CSV file.

    '''

    hits_mgr = hits.BestHits()

    df = MafParser(maf_fn).read()
    df = hits_mgr.best_hits(df)
    df.to_csv(output_fn, index=False)

