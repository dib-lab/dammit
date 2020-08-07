# Copyright (C) 2015-2019 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import logging
import os
import sys
import yaml

import click
from click import get_current_context

# TODO: put cloup on conda so we don't have to vendor it
from dammit import cloup

from dammit.log import start_logging
from dammit.config import CONFIG
from dammit.meta import __version__, __authors__, __description__, __year__

from dammit.components.convert import (maf_to_gff3_cmd,
                                       shmlast_to_gff3_cmd,
                                       hmmscan_to_gff3_cmd,
                                       cmscan_to_gff3_cmd)
from dammit.components.fastx import (rename_fasta_cmd,
                                     transcriptome_stats_cmd,
                                     annotate_fasta_cmd)
from dammit.components.filter import maf_best_hits_cmd
from dammit.components.gff3 import merge_gff3_cmd
from dammit.components.hmmer import remap_hmmer_coords_cmd
from dammit.components.run import run_group, annotate_cmd, databases_cmd
from dammit.components.config import config_group

banner = f'''
\b
     _                           _ _   
  __| | __ _ _ __ ___  _ __ ___ (_) |_ 
 / _` |/ _` | '_ ` _ \| '_ ` _ \| | __|
| (_| | (_| | | | | | | | | | | | | |_ 
 \__,_|\__,_|_| |_| |_|_| |_| |_|_|\__|
                                       
\b
{__description__}
\b
v{__version__}, {__year__}
by {" and ".join(__authors__)}
'''


@cloup.group(help=banner,
             align_sections=True)
@click.version_option(version=__version__, message='%(version)s')
@click.pass_context
def main(ctx):
    logger = logging.getLogger('dammit.component')
    start_logging()
    ctx.obj = CONFIG


main.section('Primary annotation and configuration commands',
    run_group,
    config_group
)

main.section('FASTA munging commands',
    rename_fasta_cmd,
    transcriptome_stats_cmd,
    annotate_fasta_cmd
)

main.section('Filtering commands',
    maf_best_hits_cmd
)

main.section('Conversion commands',
    maf_to_gff3_cmd,
    shmlast_to_gff3_cmd,
    hmmscan_to_gff3_cmd,
    cmscan_to_gff3_cmd
)

main.section('Transformation commands',
    merge_gff3_cmd,
    remap_hmmer_coords_cmd
)