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
