# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import csv
from itertools import count
import pandas as pd
import sys
import warnings

from dammit.fileio.base import convert_dtypes, ChunkParser, EmptyFile, warn_empty
from dammit.utils import touch

gff_version = '3.2.1'

def id_gen_wrapper():
    ID_GEN = count()
    def get_next():
        return next(ID_GEN)
    return get_next

next_ID = id_gen_wrapper()

class GFF3Parser(ChunkParser):

    columns = [('seqid', str),
                ('source', str),
                ('type', str),
                ('start', int),
                ('end', int),
                ('score', float),
                ('strand', str),
                ('phase', float),
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

    def empty(self):
        df = super(GFF3Parser, self).empty()
        df['ID'] = None
        df['Name'] = None
        df['Target'] = None
        return df

    def __iter__(self):
        '''Yields DataFrames of length chunksize from a given
        GTF/GFF file.

        GFF3 uses a 1-based, fully closed interval. Truly the devil's format.

        We convert to proper 0-based, half-open intervals.

        Yields:
            DataFrame: Pandas DataFrame with the results.
        '''
        # Read everything into a DataFrame
        n_chunks = 0
        for group in  pd.read_table(self.filename, delimiter='\t', comment='#',
                                    names=[k for k,_ in self.columns], na_values='.',
                                    converters={'attributes': self.decompose_attr_column},
                                    chunksize=self.chunksize, header=None,
                                    dtype=dict(self.columns)):

            # Generate a new DataFrame from the attributes dicts, and merge it in
            group.reset_index(drop=True, inplace=True)
            df = pd.merge(group, pd.DataFrame(list(group.attributes)),
                          left_index=True, right_index=True)
            del df['attributes']

            # Repent, repent!
            df.start = df.start - 1

            n_chunks += 1
            yield df


def maf_to_gff3(maf_df, tag='', database='',
                   ftype='translated_nucleotide_match'):
    '''Convert a MAF DataFrame to a GFF3 DataFrame ready to be written to disk.

    Args:
        maf_df (pandas.DataFrame): The MAF DataFrame. See
            dammit.fileio.maf.MafParser for column specs.
        tag (str): Extra tag to add to the source column.
        database (str): For the database entry in the attributes column.
        ftype (str): The feature type; GMOD compliant if possible.
    
    Returns:
        pandas.DataFrame: The GFF3 compliant DataFrame.
    '''

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = maf_df['q_name']
    source = '{0}.LAST'.format(tag) if tag else 'LAST'
    gff3_df['source'] = [source] * len(maf_df)
    gff3_df['type'] = [ftype] * len(maf_df)
    gff3_df['start'] = maf_df['q_start']
    gff3_df['end'] = maf_df['q_start'] + maf_df['q_aln_len']
    gff3_df['score'] = maf_df['E']
    gff3_df['strand'] = maf_df['q_strand']
    gff3_df['phase'] = ['.'] * len(maf_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next_ID()))
        data.append('Name={0}'.format(row.s_name))
        data.append('Target={0} {1} {2} {3}'.format(row.s_name, row.s_start,
                                                 row.s_start + row.s_aln_len,
                                                 row.s_strand))
        if database:
            data.append('database={0}'.format(database))

        return ';'.join(data)

    gff3_df['attributes'] = maf_df.apply(build_attr, axis=1)
    GFF3Writer.mangle_coordinates(gff3_df)

    return gff3_df


def shmlast_to_gff3(df, database=''):
    return maf_to_gff3(df, tag='shmlast', database=database,
                          ftype='conditional_reciprocal_best_LAST')


def hmmscan_to_gff3(hmmscan_df, tag='', database=''):

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = hmmscan_df['query_name']
    source = '{0}.HMMER'.format(tag) if tag else 'HMMER'
    gff3_df['source'] = [source] * len(hmmscan_df)
    gff3_df['type'] = ['protein_hmm_match'] * len(hmmscan_df)

    gff3_df['start'] = hmmscan_df['ali_coord_from']
    gff3_df['end'] = hmmscan_df['ali_coord_to']

    # Confirm whether this is the appropriate value to use
    gff3_df['score'] = hmmscan_df['domain_i_evalue']
    gff3_df['strand'] = ['.'] * len(hmmscan_df)
    gff3_df['phase'] = ['.'] * len(hmmscan_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next_ID()))
        data.append('Name={0}'.format(row.target_name))
        data.append('Target={0} {1} {2} +'.format(row.target_name,
                                                    row.hmm_coord_from+1,
                                                    row.hmm_coord_to))
        data.append('Note={0}'.format(row.description))
        data.append('accuracy={0}'.format(row.accuracy))
        data.append('env_coords={0} {1}'.format(row.env_coord_from+1,
                                                row.env_coord_to))
        if database:
            data.append('Dbxref="{0}:{1}"'.format(database,
                                                  row.target_accession))
        return ';'.join(data)

    gff3_df['attributes'] = hmmscan_df.apply(build_attr, axis=1)
    GFF3Writer.mangle_coordinates(gff3_df)

    return gff3_df


def cmscan_to_gff3(cmscan_df, tag='', database=''):

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = cmscan_df['query_name']
    source = '{0}.Infernal'.format(tag) if tag else 'Infernal'
    gff3_df['source'] = [source] * len(cmscan_df)

    # For now, using:
    #  http://www.sequenceontology.org/browser/current_svn/term/SO:0000122
    # There are more specific features for secondary structure which should
    # be extracted from the Rfam results eventually
    gff3_df['type'] = ['RNA_sequence_secondary_structure'] * len(cmscan_df)

    gff3_df['start'] = cmscan_df['seq_from']
    gff3_df['end'] = cmscan_df['seq_to']
    gff3_df['score'] = cmscan_df['e_value']
    gff3_df['strand'] = cmscan_df['strand']
    gff3_df['phase'] = ['.'] * len(cmscan_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next_ID()))
        data.append('Name={0}'.format(row.target_name))
        data.append('Target={0} {1} {2} +'.format(row.target_name,
                                                    row.mdl_from+1,
                                                    row.mdl_to))
        data.append('Note={0}'.format(row.description))
        if database:
            data.append('Dbxref="{0}:{1}"'.format(database,
                                                  row.target_accession))
        data.append('trunc={0}'.format(row.trunc))
        data.append('bitscore={0}'.format(row.score))

        return ';'.join(data)

    gff3_df['attributes'] = cmscan_df.apply(build_attr, axis=1)
    GFF3Writer.mangle_coordinates(gff3_df)

    return gff3_df


class GFF3Writer(object):

    version_line = '##gff-version 3.2.1'

    def __init__(self, filename=None, converter=None, **converter_kwds):
        self.filename = filename
        self.converter = converter
        self.converter_kwds = converter_kwds
        self.created = False

    def convert(self, data_df):
        return self.converter(data_df, **self.converter_kwds)

    def write(self, data_df, version_line=True):
        '''Write the given data to a GFF3 file, using the converter if given.

        Generates an empty file if given an empty DataFrame.

        Args:
            version_line (bool): If True, write the GFF3 version line at the.
            Note that this will cause an existing file to be overwritten, but
            will only be added in the first call to `write`.
        '''

        if self.filename is None:
            raise ValueError('Trying to write to filename None! Give GFF3Writer'
                             ' a filename.')

        if len(data_df) == 0:
            warn_empty('Writing out an empty GFF3 file to {0}'.format(self.filename))
            touch(self.filename)
            return
            
        if not self.created and version_line is True:
            with open(self.filename, 'w') as fp:
                fp.write(self.version_line + '\n')

        with open(self.filename, 'a') as fp:
            self.created = True

            if self.converter is not None:
                data_df = self.convert(data_df)
            else:
                self.mangle_coordinates(data_df)

            data_df.to_csv(fp, sep='\t', na_rep='.', columns=[k for k, v in GFF3Parser.columns],
                           index=False, header=False, quoting=csv.QUOTE_NONE,
                           float_format='%.6e')

    @staticmethod
    def mangle_coordinates(gff3_df):
        '''Although 1-based fully closed intervals are of the Beast,
        we will respect the convention in the interests of peace between
        worlds and compatibility.

        Args:
            gff3_df (DataFrame): The DataFrame to "fix".
        '''
        gff3_df.start += 1


