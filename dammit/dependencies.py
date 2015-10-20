#!/usr/bin/env python
from __future__ import print_function

import logging
import os
import sys

from doit.dependency import Dependency, SqliteDB

from common import which
import common
from tasks import get_download_and_untar_task

logger = logging.getLogger(__name__)


def check_system_path():
    
    deps = {}

    hmmscan = which('hmmscan')
    hmmpress = which('hmmpress')
    if hmmscan is None or hmmpress is None:
        deps['HMMER'] = False
    else:
        deps['HMMER'] = True
        logger.debug('hmmscan:' + hmmscan)
        logger.debug('hmmpress:' + hmmpress)

    cmscan = which('cmscan')
    cmpress = which('cmpress')
    if cmscan is None or cmpress is None:
        deps['Infernal'] = False
    else:
        deps['Infernal'] = True
        logger.debug('cmscan:' + cmscan)
        logger.debug('cmpress:' + cmpress)

    blastp = which('blastp')
    blastx = which('blastx')
    tblastn = which('tblastn')
    makeblastdb = which('makeblastdb')
    if (blastp is None) or (blastx is None) \
        or (tblastn is None) or (makeblastdb is None):
        
        deps['BLAST+'] = False
    else:
        deps['BLAST+'] = True
        logger.debug('blastp:' + blastp)
        logger.debug('blastx:' + blastx)
        logger.debug('tblastn:' + tblastn)
        logger.debug('makeblastdb:' + makeblastdb)

    busco = which('BUSCO_v1.1b1.py')
    if busco is None:
        deps['BUSCO'] = False
    else:
        deps['BUSCO'] = True
        logger.debug('BUSCO:' + busco)


    longorfs = which('TransDecoder.LongOrfs')
    predict = which('TransDecoder.Predict')
    if longorfs is None or predict is None:
        deps['TransDecoder'] = False
    else:
        deps['TransDecoder'] = True
        logger.debug('TransDecoder.LongOrfs:' + longorfs)
        logger.debug('TransDecoder.Predict:' + predict)

    lastdb = which('lastdb')
    lastal = which('lastal')
    if lastdb is None or lastal is None:
        deps['LAST'] = False
    else:
        deps['LAST'] = True
        logger.debug('lastal:' + lastal)
        logger.debug('latsdb:' + lastdb)

    crb_blast = which('crb-blast')
    if crb_blast is None:
        deps['crb-blast'] = False
    else:
        deps['crb-blast'] = True
        logger.debug('crb-blast:' + crb_blast)

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

