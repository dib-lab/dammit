#!/usr/bin/env python
from __future__ import print_function

from .common import which
from .tasks import doit_task

import os


def parallel_fasta(input_filename, n_jobs, file_size=10000):
    exc = which('parallel')
    block_size = int(file_size) / int(n_jobs)
    cmd = ['cat', input_filename, '|', exc, '--block', str(block_size),
           '--pipe', '--recstart', '">"', '--gnu', '-j', str(n_jobs)]

    return ' '.join(cmd)

def multinode_parallel_fasta(input_filename, ppn, nodes, file_size):
    '''
    exports = 'export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH '\
              '--env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION '\
              '--env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME '\
              '--env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED '\
              '--env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE '\
              '--env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY";'
    '''
    exc = which('parallel')
    block_size = int(file_size) / int(ppn * nodes)
    cmd = ['cat', input_filename, '|', exc, '--block', str(block_size),
           '--pipe', '--recstart', '">"', '--gnu', '--jobs 1', 
           '--sshloginfile $PBS_NODEFILE', '--workdir $PWD']

    return ' '.join(cmd)


@doit_task
def get_filesize_task(filename):

    def get_filesize(filename):
        return {'size': os.path.getsize(filename)}

    return {'name': 'get_filesize:{0}'.format(filename),
            'actions': [(get_filesize, [filename])],
            'file_dep': [filename]}

