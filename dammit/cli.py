# Copyright (C) 2015-2019 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

from copy import deepcopy
import os

import click

# TODO: put cloup on conda so we don't have to vendor it
from dammit import cloup

from dammit.config import get_config_obj
from dammit.meta import __version__, __authors__, __description__, __year__, __path__
from dammit.utils import Namespace

from dammit.components.convert import (maf_to_gff3_cmd,
                                       shmlast_to_gff3_cmd,
                                       hmmscan_to_gff3_cmd,
                                       cmscan_to_gff3_cmd,
                                       busco_to_gff3_cmd)
from dammit.components.fastx import (rename_fasta_cmd,
                                     transcriptome_stats_cmd,
                                     annotate_fasta_cmd)
from dammit.components.filter import maf_best_hits_cmd
from dammit.components.gff3 import merge_gff3_cmd
from dammit.components.hmmer import remap_hmmer_coords_cmd
from dammit.components.run import run_group
from dammit.components.config import config_group

banner_txt = open(os.path.join(__path__, 'banner.txt')).read().rstrip('\n')

banner = f'''
\b
{banner_txt}                                      
{__description__}
\b
v{__version__}, {__year__}
by {" and ".join(__authors__)}
'''


@cloup.group(help=banner,
             align_sections=True)
@click.version_option(version=__version__, message='%(version)s')
@click.pass_context
@click.option('--config-file',
              help='A YAML or JSON file providing values to override'\
                   ' built-in config. Advanced use only!')
def main(ctx, config_file):
    CONFIG = get_config_obj(config_file)
    CONFIG.banner = banner
    CONFIG.gui = Namespace()
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
    cmscan_to_gff3_cmd,
    busco_to_gff3_cmd
)


main.section('Transformation commands',
    merge_gff3_cmd,
    remap_hmmer_coords_cmd
)