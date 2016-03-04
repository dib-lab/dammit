#!/usr/bin/env python

import csv
from itertools import count
import pandas as pd

gff3_cols = [('seqid', str),
             ('source', str),
             ('type', str),
             ('start', int),
             ('end', int),
             ('score', float),
             ('strand', str),
             ('phase', str),
             ('attributes', str)]

gff_version = '3.2.1'

ID_GEN = count()

def version_line():
    return '##gff-version {v}'.format(v=gff_version)


def write_gff3_df(df, fp):

    df.to_csv(fp, sep='\t', na_rep='.', columns=[k for k, v in gff3_cols],
              index=False, header=False, quoting=csv.QUOTE_NONE)

def mangle_coordinates(gff3_df):
    '''Although 1-based fully closed intervals are of the Beast,
    we will respect the convention in the interest of peace between
    worlds and compatibility.

    Args:
        gff3_df (DataFrame): The DataFrame to "fix".
    '''
    gff3_df.start += 1

def maf_to_gff3_df(maf_df, tag='', database=''):

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = maf_df['q_name']
    source = '{0}.LAST'.format(tag) if tag else 'LAST'
    gff3_df['source'] = [source] * len(maf_df)
    gff3_df['type'] = ['translated_nucleotide_match'] * len(maf_df)
    gff3_df['start'] = maf_df['q_start']
    gff3_df['end'] = maf_df['q_start'] + maf_df['q_aln_len']
    gff3_df['score'] = maf_df['E']
    gff3_df['strand'] = maf_df['q_strand']
    gff3_df['phase'] = ['.'] * len(maf_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
        data.append('Name={0}'.format(row.s_name))
        data.append('Target={0} {1} {2} {3}'.format(row.s_name, row.s_start,
                                                 row.s_start + row.s_aln_len,
                                                 row.s_strand))
        if database:
            data.append('database={0}'.format(database))

        return ';'.join(data)

    gff3_df['attributes'] = maf_df.apply(build_attr, axis=1)
    mangle_coordinates(gff3_df)

    return gff3_df


def crb_to_gff3_df(crb_df, tag='', database=''):

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = crb_df['query']
    source = '{0}.crb-blast'.format(tag) if tag else 'crb-blast'
    gff3_df['source'] = [source] * len(crb_df)
    gff3_df['type'] = ['translated_nucleotide_match'] * len(crb_df)
    gff3_df['start'] = crb_df['qstart']
    gff3_df['end'] = crb_df['qend']
    gff3_df['score'] = crb_df['evalue']
    gff3_df['strand'] = crb_df['qstrand']
    gff3_df['phase'] = ['.'] * len(crb_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
        data.append('Name={0}'.format(row.subject))
        data.append('Target={0} {1} {2} {3}'.format(row.subject, row.sstart,
                                                    row.send, row.sstrand))

        if database:
            data.append('database={0}'.format(database))

        return ';'.join(data)

    gff3_df['attributes'] = crb_df.apply(build_attr, axis=1)
    mangle_coordinates(gff3_df)

    return gff3_df


def hmmscan_to_gff3_df(hmmscan_df, tag='', database=''):

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
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
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
    mangle_coordinates(gff3_df)

    return gff3_df


def cmscan_to_gff3_df(cmscan_df, tag='', database=''):

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
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
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
    mangle_coordinates(gff3_df)

    return gff3_df
