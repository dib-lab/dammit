# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

from .utils import run, run_shell_cmd


def test_dammit_version():
    '''Test the dammit --version command.
    '''

    from dammit.meta import __version__
    status, out, err = run('--version')
    print(status, out, err)
    assert status == 0
    assert out.strip() ==  __version__


class TestSnakemakeExtraArgs:
    
    def test_annotate_cmd(self, tmpdir, datadir):
        '''Test passing extra args to be consumed by Snakemake.
        '''
        
        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')

            args = ['dammit', 'run', 
                    '--busco-group', 'bacteria_odb10',
                    '--busco-group', 'saccharomycetes_odb10',
                    '--pipeline', 'quick',
                    'annotate', transcripts, '--dag']
            status, out, err = run_shell_cmd(' '.join(args))
            
            assert status == 0
            assert 'Beginning workflow execution...' in err
            assert out.startswith('digraph snakemake_dag {')