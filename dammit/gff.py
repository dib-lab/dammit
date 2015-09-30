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

ID_GEN = count()

def write_gff3_df(df, fp):

    df.to_csv(fp, sep='\t', na_rep='.', columns=[k for k, v in gff3_cols],
              index=False, header=False, quoting=csv.QUOTE_NONE)

def maf_to_gff3_df(maf_df, tag, database=''):

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = maf_df['q_name']
    gff3_df['source'] = ['LAST'] * len(maf_df)
    gff3_df['type'] = ['translated_nucleotide_match'] * len(maf_df)
    gff3_df['start'] = maf_df['q_start'] + 1
    gff3_df['end'] = maf_df['q_start'] + maf_df['q_aln_len'] + 1
    gff3_df['score'] = maf_df['E']
    gff3_df['strand'] = maf_df['q_strand']
    gff3_df['phase'] = ['.'] * len(maf_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
        data.append('Name={0}:{1}'.format(row.s_name, tag))
        data.append('Target={0} {1} {2} {3}'.format(row.s_name, row.s_start,
                                                 row.s_start + row.s_aln_len,
                                                 row.s_strand))
        if database:
            data.append('database={0}'.format(database))

        return ';'.join(data)

    gff3_df['attributes'] = maf_df.apply(build_attr, axis=1)
    
    return gff3_df

def blast_to_gff3_df(blast_df, prog, RBH=False, database=''):

    assert prog in ['BLASTX', 'TBLASTN', 'BLASTP']

    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = blast_df['qseqid']
    gff3_df['source'] = [prog] * len(blast_df)

    ftype = 'protein_match'
    if prog in ['BLASTX', 'TBLASTN']:
        ftype = 'translated_nucleotide_match'
    gff3_df['type'] = [ftype] * len(blast_df)

    gff3_df['start'] = blast_df['qstart'] + 1
    gff3_df['end'] = blast_df['qend']
    gff3_df['score'] = blast_df['evalue']
    gff3_df['strand'] = blast_df['qstrand']
    gff3_df['phase'] = ['.'] * len(blast_df)

    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
        data.append('Name={0}:{1}'.format(row.sseqid, tag))
        data.append('Target={0} {1} {2} {3}'.format(row.sseqid, row.sstart,
                                                    row.send, row.sstrand))
        if database:
            data.append('database={0}'.format(database))

        return ';'.join(data)

    gff3_df['attributes'] = blast_df.apply(build_attr, axis=1)
    return gff3_df

def hmmscan_to_gff3_df(hmmscan_df, tag, database=''):
    
    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = hmmscan_df['query_name']
    gff3_df['source'] = ['HMMER'] * len(hmmscan_df)
    gff3_df['type'] = ['protein_hmm_match'] * len(hmmscan_df)

    # This is kludgy and quite wrong
    gff3_df['start'] = hmmscan_df['env_coord_from']
    gff3_df['end'] = hmmscan_df['env_coord_to'] * 3

    # Confirm whether this is the appropriate value to use
    gff3_df['score'] = hmmscan_df['domain_i_evalue']
    gff3_df['strand'] = ['.'] * len(hmmscan_df)
    gff3_df['phase'] = ['.'] * len(hmmscan_df)
    
    def build_attr(row):
        data = []
        data.append('ID=homology:{0}'.format(next(ID_GEN)))
        data.append('Name={0}:{1}'.format(row.target_name, tag))
        data.append('Target={0} {1} {2} +'.format(row.target_name,
                                                    row.hmm_coord_from,
                                                    row.hmm_coord_to))
        data.append('Note={0}'.format(row.description))
        data.append('accuracy={0}'.format(row.accuracy))
        if database:
            data.append('Dbxref="{0}:{1}"'.format(database,
                                                  row.target_accession))
        return ';'.join(data)

    gff3_df['attributes'] = hmmscan_df.apply(build_attr, axis=1)

    return gff3_df

def cmscan_to_gff3_df(cmscan_df, tag, database=''):
    
    gff3_df = pd.DataFrame()
    gff3_df['seqid'] = cmscan_df['query_name']
    gff3_df['source'] = ['Infernal'] * len(cmscan_df)

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
        data.append('Name={0}:{1}'.format(row.target_name, tag))
        data.append('Target={0} {1} {2} +'.format(row.target_name,
                                                    row.mdl_from,
                                                    row.mdl_to))
        data.append('Note={0}'.format(row.description))
        if database:
            data.append('Dbxref="{0}:{1}"'.format(database,
                                                  row.target_accession))
        data.append('trunc={0}'.format(row.trunc))
        data.append('bitscore={0}'.format(row.score))

        return ';'.join(data)

    gff3_df['attributes'] = cmscan_df.apply(build_attr, axis=1)

    return gff3_df

