# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import pandas as pd
from dammit.fileio.base import convert_dtypes, ChunkParser


class InfernalParser(ChunkParser):

    columns = [('target_name', str),
               ('target_accession', str),
               ('query_name', str),
               ('query_accession', str),
               ('mdl', str),
               ('mdl_from', int),
               ('mdl_to', int),
               ('seq_from', int),
               ('seq_to', int),
               ('strand', str),
               ('trunc', str),
               ('pass', str),
               ('gc', float),
               ('bias', float),
               ('score', float),
               ('e_value', float),
               ('inc', str),
               ('description', str)]
    
    def __init__(self, filename, **kwargs):
        super(InfernalParser, self).__init__(filename, **kwargs)

    def __iter__(self):
        '''Yields DataFrames of length chunksize from a given
        cmscan result file.

        The format uses 1-based, fully open intervals; when the strand
        is negative, the start coordinate is larger than the end. 
        Truly Infernal.

        We convert to proper 0-based, half-open, ordered intervals.

        Yields:
            DataFrame: Pandas DataFrame with the cmscan hits.
        '''

        data = []
        n_entries = 0
        with open(self.filename) as fp:
            for ln in fp:
                ln = ln.strip()
                if not ln or ln.startswith('#'):
                    continue
                tokens = ln.split()
                data.append(tokens[:len(self.columns)-1] + \
                            [' '.join(tokens[len(self.columns)-1:])])
                n_entries += 1
                if len(data) >= self.chunksize:
                    yield self._build_df(data)
                    data = []

        if n_entries == 0:
            self.raise_empty()
        if data:
            yield self._build_df(data)

    def _build_df(self, data):
        if not data:
            self.raise_empty()
        df = pd.DataFrame(data, columns=[k for k, _ in self.columns])
        convert_dtypes(df, dict(self.columns))
        # fix the evil coordinate system
        sidx = df.seq_from > df.seq_to
        df.loc[sidx, 'seq_from'], df.loc[sidx, 'seq_to'] = \
            df.loc[sidx, 'seq_to'], df.loc[sidx, 'seq_from']
        df.mdl_from = df.mdl_from - 1
        df.seq_from = df.seq_from - 1
        return df
