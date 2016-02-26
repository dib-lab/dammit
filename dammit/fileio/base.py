class BaseParser(object):

    def __init__(self, filename):
        self.filename = filename

class ChunkParser(BaseParser):

    def __init__(self, filename, chunksize=10000):
        super(BaseParser, self).__init__(filename)
        self.chunksize = chunksize

    def __iter__(self):
        raise NotImplementedError()
        yield
