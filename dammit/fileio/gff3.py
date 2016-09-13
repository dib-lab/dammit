import pandas as pd
from .base import convert_dtypes, ChunkParser

class GFF3Parser(ChunkParser):

    columns = [('seqid', str),
                ('source', str),
                ('feature_type', str),
                ('start', float),
                ('end', float),
                ('score', float),
                ('strand', str),
                ('frame', float),
                ('attributes', str)]
    
    def __init__(self, filename, **kwargs):
        super(GFF3Parser, self).__init__(filename, **kwargs)

    @staticmethod
    def decompose_attr_column(col):
        d = {}
        for item in col.strip(';').split(';'):
            key, _, val = item.strip().partition('=')
            d[key] = val.strip('')
        return d

    def __iter__(self):
        '''Yields DataFrames of length chunksize from a given
        GTF/GFF file.

        GFF3 uses a 1-based, fully closed interval. Truly the devil's format.

        We convert to proper 0-based, half-open intervals.

        Yields:
            DataFrame: Pandas DataFrame with the results.
        '''
        # Read everything into a DataFrame
        for group in  pd.read_table(self.filename, delimiter='\t', comment='#',
                                    names=[k for k,_ in self.columns], na_values='.',
                                    converters={'attributes': self.decompose_attr_column},
                                    chunksize=chunksize, header=None,
                                    dtype=dict(gff_cols)):

            # Generate a new DataFrame from the attributes dicts, and merge it in
            df = pd.merge(group, pd.DataFrame(list(group.attributes)),
                          left_index=True, right_index=True)
            del df['attributes']

            # Repent, repent!
            df.start = df.start - 1

            yield df
