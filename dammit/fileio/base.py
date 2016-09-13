from itertools import count


def next_or_raise(fp):
    counter = count()
    def func(raise_exc=True):
        line = fp.readline()
        n = next(counter)
        if raise_exc is True and line == '':
            raise RuntimeError('Malformed MAF file (line {0})'.format(n))
        return line
    return func


def convert_dtypes(df, dtypes):
    '''Convert the columns of a DataFrame to the types specified
    in the given dictionary, inplace.

    Args:
        df (DataFrame): The DataFrame to convert.
        dtypes (dict): Dictionary mapping columns to types.
    '''

    for c in df.columns:
        try:
            df[c] = df[c].astype(dtypes[c])
        except KeyError:
            pass


class BaseParser(object):

    def __init__(self, filename):
        self.filename = filename


class ChunkParser(BaseParser):

    def __init__(self, filename, chunksize=10000):
        self.chunksize = chunksize
        super(ChunkParser, self).__init__(filename)

    def __iter__(self):
        raise NotImplementedError()
        yield

    def read(self):
        '''Read the entire file at once and return as a single DataFrame.
        '''
        return pd.concat(self, ignore_index=True)

    def empty(self):
        df = pd.DataFrame(columns=[k for k, _ in self.columns])
        convert_dtypes(df, dict(self.columns))
        return df
