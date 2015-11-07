#!/usr/bin/env python
from __future__ import print_function

import logging
import os
import sys

from doit.dependency import Dependency, SqliteDB

from .common import which
from . import common
from .tasks import get_download_and_untar_task


class DependencyHandler(object):

    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)

    def handle(self):

        common.print_header('Checking PATH for dependencies', level=2)

        system_deps = self.check_system_path()
        
        missing = []
        for key, status in system_deps.iteritems():
            if status is False:
                missing.append(key)
                self.logger.warning('[ ] {0}'.format(key))
            else:
                self.logger.info('[x] {0}'.format(key))

        common.print_header('Dependency results', level=2)

        if missing:
            self.logger.warning('{0} missing'.format(', '.join(missing)))
        else:
            self.logger.info('All dependencies satisfied!')

        return missing

    def check_or_fail(self):
        missing = self.handle()
        if missing:
            self.logger.error('Install dependencies to continue; exiting')
            sys.exit(1)

    def check_system_path(self):
        
        deps = {}

        hmmscan = which('hmmscan')
        hmmpress = which('hmmpress')
        if hmmscan is None or hmmpress is None:
            deps['HMMER'] = False
        else:
            deps['HMMER'] = True
            self.logger.debug('hmmscan:' + hmmscan)
            self.logger.debug('hmmpress:' + hmmpress)

        cmscan = which('cmscan')
        cmpress = which('cmpress')
        if cmscan is None or cmpress is None:
            deps['Infernal'] = False
        else:
            deps['Infernal'] = True
            self.logger.debug('cmscan:' + cmscan)
            self.logger.debug('cmpress:' + cmpress)

        blastp = which('blastp')
        blastx = which('blastx')
        tblastn = which('tblastn')
        makeblastdb = which('makeblastdb')
        if (blastp is None) or (blastx is None) \
            or (tblastn is None) or (makeblastdb is None):
            
            deps['BLAST+'] = False
        else:
            deps['BLAST+'] = True
            self.logger.debug('blastp:' + blastp)
            self.logger.debug('blastx:' + blastx)
            self.logger.debug('tblastn:' + tblastn)
            self.logger.debug('makeblastdb:' + makeblastdb)

        busco = which('BUSCO_v1.1b1.py')
        if busco is None:
            deps['BUSCO'] = False
        else:
            deps['BUSCO'] = True
            self.logger.debug('BUSCO:' + busco)


        longorfs = which('TransDecoder.LongOrfs')
        predict = which('TransDecoder.Predict')
        if longorfs is None or predict is None:
            deps['TransDecoder'] = False
        else:
            deps['TransDecoder'] = True
            self.logger.debug('TransDecoder.LongOrfs:' + longorfs)
            self.logger.debug('TransDecoder.Predict:' + predict)

        lastdb = which('lastdb')
        lastal = which('lastal')
        if lastdb is None or lastal is None:
            deps['LAST'] = False
        else:
            deps['LAST'] = True
            self.logger.debug('lastal:' + lastal)
            self.logger.debug('latsdb:' + lastdb)

        crb_blast = which('crb-blast')
        if crb_blast is None:
            deps['crb-blast'] = False
        else:
            deps['crb-blast'] = True
            self.logger.debug('crb-blast:' + crb_blast)

        return deps


