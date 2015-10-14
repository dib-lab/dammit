#!/usr/bin/env python
from __future__ import print_function

import logging
import os
from platform import system
import sys

from doit.dependency import Dependency, SqliteDB

import common
from tasks import get_download_and_untar_task

logger = logging.getLogger(__name__)

def which(program):
    '''Checks whether the given program (or program path) is valid and
    executable.

    NOTE: Sometimes copypasta is okay! This function came from stackoverflow:

        http://stackoverflow.com/a/377028/5109965

    Args:
        program (str): Either a program name or full path to a program.

    Returns:
        Return the path to the executable or None if not found
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def check_system_path():
    
    deps = {}

    hmmscan = which('hmmscan')
    hmmpress = which('hmmpress')
    logger.debug('hmmscan:' + hmmscan)
    logger.debug('hmmpress:' + hmmpress)
    if hmmscan is None or hmmpress is None:
        deps['HMMER'] = False
    else:
        deps['HMMER'] = True

    cmscan = which('cmscan')
    cmpress = which('cmpress')
    logger.debug('cmscan:' + cmscan)
    logger.debug('cmpress:' + cmpress)
    if cmscan is None or cmpress is None:
        deps['Infernal'] = False
    else:
        deps['Infernal'] = True

    blastp = which('blastp')
    blastx = which('blastx')
    tblastn = which('tblastn')
    makeblastdb = which('makeblastdb')
    logger.debug('blastp:' + blastp)
    logger.debug('blastx:' + blastx)
    logger.debug('tblastn:' + tblastn)
    logger.debug('makeblastdb:' + makeblastdb)
    if (blastp is None) or (blastx is None) \
        or (tblastn is None) or (makeblastdb is None):

        deps['BLAST+'] = False
    else:
        deps['BLAST+'] = True

    busco = which('BUSCO_v1.1b1.py')
    logger.debug('BUSCO:' + busco)
    if busco is None:
        deps['BUSCO'] = False
    else:
        deps['BUSCO'] = True


    longorfs = which('TransDecoder.LongOrfs')
    predict = which('TransDecoder.Predict')
    logger.debug('TransDecoder.LongOrfs:' + longorfs)
    logger.debug('TransDecoder.Predict:' + predict)
    if longorfs is None or predict is None:
        deps['TransDecoder'] = False
    else:
        deps['TransDecoder'] = True

    lastdb = which('lastdb')
    lastal = which('lastal')
    logger.debug('lastal:' + lastal)
    logger.debug('latsdb:' + lastdb)
    if lastdb is None or lastal is None:
        deps['LAST'] = False
    else:
        deps['LAST'] = True

    crb_blast = which('crb-blast')
    logger.debug('crb-blast:' + crb_blast)
    if crb_blast is None:
        deps['crb-blast'] = False
    else:
        deps['crb-blast'] = True

    return deps


def do_check():

    common.print_header('Checking PATH for dependencies', level=2)

    system_deps = check_system_path()
    
    missing = []
    for key, status in system_deps.iteritems():
        if status is False:
            missing.append(key)
            logger.warning('[ ] {0}'.format(key))
        else:
            logger.info('[x] {0}'.format(key))

    common.print_header('Dependency results', level=2)

    if missing:
        logger.warning('{0} missing'.format(', '.join(missing)))
    else:
        logger.info('All dependencies satisfied!')

    return missing

