#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : run.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.10.2020

from ope.io.gff3 import GFF3Converter


class MAF_to_GFF3(GFF3Converter):

    def __init__(self, maf_df, tag='', database='',
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
        super().__init__(maf_df, tag, database, ftype)

    def seqid(self):
        return self.from_df['q_name']

    def source(self):
        source = '{0}.LAST'.format(self.tag) if self.tag else 'LAST'
        return [source] * len(self.from_df)
    
    def feature_type(self):
        return [self.ftype] * len(self.from_df)

    def start(self):
        return self.from_df['q_start']
    
    def end(self):
        return self.from_df['q_start'] + self.from_df['q_aln_len']
    
    def score(self):
        return self.from_df['E']
    
    def strand(self):
        return self.from_df['q_strand']
    
    def phase(self):
        return ['.'] * len(self.from_df)
    
    def ID_attr(self, IDs):
        return 'homology:' + IDs

    def attr_from_row(self, row):
        attrs = {'Name': '{0}'.format(row.s_name),
                 'Target': '{0} {1} {2} {3}'.format(row.s_name, row.s_start + 1,
                                                    row.s_start + row.s_aln_len,
                                                    row.s_strand)}
        if self.database:
            attrs['database'] = '{0}'.format(self.database)

        return attrs


class Shmlast_to_GFF3(MAF_to_GFF3):

    def __init__(self, shmlast_df, tag='shmlast', database='',
                       ftype='conditional_reciprocal_best_LAST'):
        super().__init__(shmlast_df, tag=tag, database=database, ftype=ftype)


class HMMScan_to_GFF3(GFF3Converter):

    def __init__(self, hmmscan_df, tag='', database=''):
        super().__init__(hmmscan_df, tag=tag, database=database,
                         ftype='protein_hmm_match')

    def seqid(self):
        return self.from_df['query_name']
    
    def source(self):
        source = '{0}.HMMER'.format(self.tag) if self.tag else 'HMMER'
        return [source] * len(self.from_df)
    
    def feature_type(self):
        return [self.ftype] * len(self.from_df)

    def start(self):
        return self.from_df['ali_coord_from']
    
    def end(self):
        return self.from_df['ali_coord_to']

    def score(self):
        # Confirm whether this is the appropriate value to use
        return self.from_df['domain_i_evalue']
    
    def strand(self):
        return ['.'] * len(self.from_df)
    
    def phase(self):
        return ['.'] * len(self.from_df)

    def ID_attr(self, IDs):
        return 'homology:' + IDs

    def attr_from_row(self, row):
        attrs = {'Name': '{0}'.format(row.target_name),
                 'Target': '{0} {1} {2} +'.format(row.target_name,
                                                  row.hmm_coord_from + 1,
                                                  row.hmm_coord_to),
                 'accuracy': '{0}'.format(row.accuracy),
                 'env_coords': '{0} {1}'.format(row.env_coord_from + 1,
                                                row.env_coord_to)}

        if self.database and (row.target_accession != '-'):
            attrs['Dbxref'] = '"{0}:{1}"'.format(self.database,
                                               row.target_accession)

        if row.description != '-':
            attrs['Note'] = '{0}'.format(row.description)
        
        return attrs


class CMScan_to_GFF3(GFF3Converter):

    def __init__(self, cmscan_df, tag='', database=''):
        super().__init__(cmscan_df, tag=tag, database=database,
                         ftype='RNA_sequence_secondary_structure')

    def seqid(self):
        return self.from_df['query_name']
    
    def source(self):
        source = '{0}.Infernal'.format(self.tag) if self.tag else 'Infernal'
        return [source] * len(self.from_df)

    def feature_type(self):
        # For now, using:
        #  http://www.sequenceontology.org/browser/current_svn/term/SO:0000122
        # There are more specific features for secondary structure which should
        # be extracted from the Rfam results eventually
        return['RNA_sequence_secondary_structure'] * len(self.from_df)

    def start(self):
        return self.from_df['seq_from']
    
    def end(self):
        return self.from_df['seq_to']
    
    def score(self):
        return self.from_df['e_value']
    
    def strand(self):
        return self.from_df['strand']
    
    def phase(self):
        return ['.'] * len(self.from_df)

    def ID_attr(self, IDs):
        return 'homology:' + IDs

    def attr_from_row(self, row):
        attrs = {'Name': '{0}'.format(row.target_name),
                 'Target': '{0} {1} {2} +'.format(row.target_name,
                                                  row.mdl_from+1,
                                                  row.mdl_to),
                 'trunc': '{0}'.format(row.trunc),
                 'bitscore': '{0}'.format(row.score)}

        if self.database and (row.target_accession != '-'):
            attrs['Dbxref'] = '"{0}:{1}"'.format(self.database,
                                                 row.target_accession)
        
        if row.description != '-':
            attrs['Note'] = '{0}'.format(row.description)
            
        return attrs


class BUSCO_to_GFF3(GFF3Converter):

    def __init__(self, busco_df, tag='BUSCO', database='', busco_version='4.1.1'):
        self.busco_version = busco_version
        super().__init__(busco_df, tag=tag, database=database,
                         ftype='BUSCO_ortholog')
    
    def seqid(self):
        return self.from_df['ID']
    
    def source(self):
        return [self.tag] * len(self.from_df)
    
    def feature_type(self):
        return [self.ftype] * len(self.from_df)
    
    def start(self):
        if self.busco_version == '5.0.0':
            return self.from_df['Start']

        return [0] * len(self.from_df)
    
    def end(self):
        if self.busco_version == '5.0.0':
            return self.from_df['End']

        return self.from_df['Length_tx']
    
    def score(self):
        return self.from_df['Score']
    
    def strand(self):
        return ['.'] * len(self.from_df)
    
    def phase(self):
        return ['.'] * len(self.from_df)
    
    def ID_attr(self, IDs):
        return 'busco:' + IDs
    
    def attr_from_row(self, row):
        attrs = {'Name': '{0}'.format(row.BUSCO_id),
                 'length': '{0}'.format(int(row.Length)),
                 'status': '{0}'.format(row.Status)} 
        
        return attrs