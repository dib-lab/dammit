#!/usr/bin/env python
from __future__ import print_function

import pandas as pd

class BestHits(object):

    def __init__(self, comparison_col='E', query_name_col='q_name', 
                 subject_name_col='s_name', query_length_col='q_len',
                 subject_length_col='q_len'):
        self.comparison_col = comparison_col
        self.query_name_col = query_name_col
        self.subject_name_col = subject_name_col
        self.query_length_col = query_length_col
        self.subject_length_col = subject_length_col

    def best_hits(self, aln_df):
        '''Get the best hit for each query in the alignment DataFrame.

        Operates in-place. Sorts the hits by query name and then e-value,
        and then uses the drop_duplicates() function to remove all but the
        first (lowest e-value) hit for each query.

        Args:
            aln_df (DataFrame): The MAF alignment DataFrame.
            comp_col (str): The column name to use for sorting hits.
        '''

        return aln_df.sort_values(
                   by=[self.query_name_col,self.comparison_col]
                ).drop_duplicates(subset=self.query_name_col)

    def reciprocal_best_hits(self, aln_df_A, aln_df_B):
        '''Given to DataFrames with reciprocal MAF alignments, get the
        reciprocal best hits.

        Uses the CRBL.best_hits function to get best hits for each DataFrame,
        and then does an inner join to find the reciprocals. The DataFrames
        *will* be modified in-place by the best_hits function!

        Args:
            aln_df_A (DataFrame): The query hits.
            aln_df_B (DataFrame): The subject hits.
            comp_col (str): The comparison column to use for best hits.
        Returns:
            DataFrame with the reciprocal best hits.
        '''

        aln_df_A = CRBL.best_hits(aln_df_A)
        aln_df_B = CRBL.best_hits(aln_df_B)

        # Join between subject A and query B
        rbh_df = pd.merge(aln_df_A, aln_df_B, how='inner', 
                          left_on=self.subject_name_col, 
                          right_on=self.query_name_col)

        # Renamed columns after join
        left_q_col = self.query_name_col + '_x'
        left_s_col = self.subject_name_col + '_x'
        left_len_col = self.query_length_col + '_x'
        left_comp_col = self.comparison_col + '_x'
        right_q_col = self.query_name_col + '_y'
        right_s_col = self.subject_name_col + '_y'
        right_len_col = self.query_len_col + '_y'

        # Select those where query A is the same as subject B
        rbh_df = rbh_df[(rbh_df[left_q_col] == rbh_df[right_s_col])]

        # Drop extra columns and rename for clarity
        del rbh_df[left_s_col]
        del rbh_df[right_s_col]
        rbh_df.rename(columns={left_q_col: self.query_name_col, 
                               right_q_col: self.subject_name_col,
                               left_comp_col: self.comparison_col,
                               left_len_col: self.query_length_col,
                               right_len_col: self.subject_length_col},
                      inplace=True)
        
        return rbh_df


