#/usr/bin/env python
from __future__ import print_function

import numpy as np
from sklearn.externals import joblib
from sklearn import svm
from sklearn import preprocessing


class gCRBB(object):

    def __init__(self, training_pep_fn, transcriptome_fn, database_pep_fn,
                 model_fn, classifier=svm.OneClassSVM, 
                 params={nu:0.1, kernel:'rbf', gamma:'0.1'},
                 lastal_cfg=common.CONFIG['settings']['last']['lastal'],
                 lastdb_cfg=common.CONFIG['settings']['last']['lastdb']):

        self.training_pep_fn = training_pep_fn
        self.transcriptome_fn = transcriptome_fn
        self.database_pep_fn = database_pep_fn

        self.train_x_db_fn = '{0}.x.{1}.maf'.format(self.training_pep_fn,
                                                    self.database_pep_fn)
        self.db_x_train_fn = '{0}.x.{1}.maf'.format(self.database_pep_fn,
                                                    self.training_pep_fn)
        self.transcriptome_x_db_fn = '{0}.x.{1}.maf'.format(self.transcriptome_fn,
                                                            self.database_pep_fn)

        self.model_fn = model_fn

        self.trained = False
        self.features = ['q_aln_len', 'E']
        self.classifier = classifier(**params)

    @staticmethod
    def best_hits(aln_df, comp_col='E'):
        '''Get the best hit for each query in the alignment DataFrame.

        Operates in-place. Sorts the hits by query name and then e-value,
        and then uses the drop_duplicates() function to remove all but the
        first (lowest e-value) hit for each query.

        Args:
            aln_df (DataFrame): The MAF alignment DataFrame.
            comp_col (str): The column name to use for sorting hits.
        '''

        aln_df.sort(columns=['q_name', comp_col], inplace=True)
        aln_df.drop_duplicates(subset='q_name', inplace=True)

    @staticmethod
    def reciprocal_best_hits(aln_df_A, aln_df_B, comp_col='E'):
        '''Given to DataFrames with reciprocal MAF alignments, get the
        reciprocal best hits.

        Uses the gCRBB.best_hits function to get best hits for each DataFrame,
        and then does an inner join to find the reciprocals. The DataFrames
        *will* be modified in-place by the best_hits function!

        Args:
            aln_df_A (DataFrame): The query hits.
            aln_df_B (DataFrame): The subject hits.
            comp_col (str): The comparison column to use for best hits.
        Returns:
            DataFrame with the reciprocal best hits.
        '''

        gCRBB.best_hits(aln_df_A)
        gCRBB.best_hits(alb_df_B)

        # Join between subject A and query B
        rbh_df = pd.merge(aln_df_A, aln_df_B, how='inner', 
                          left_on='s_name', right_on='q_name')
        # Select those where query A is the same as subject B
        rbh_df = rbh_df[(rbh_df['q_name_x'] == rbh_df['s_name_y'])]

        # Drop extra columns and rename for clarity
        del rbh_df['s_name_x']
        del rbh_df['s_name_y']
        rbh_df.rename(columns={'q_name_x': 'q_name', 'q_name_y': 's_name'},
                      inplace=True)
        
        return rbh_df

    def format_training_pep_task(self):
        pass

    def format_database_pep_task(self):
        pass

    def align_training_pep_task(self):
        pass

    def align_database_pep_task(self):
        pass

    def align_transcriptome_task(self):
        pass

    def preprocess(self, alignments):

        # Subset out the features
        data = alignments.dropna()[self.features].astype(float)

        # e-values of 0.0 mess things up, set them to something really low
        data['E'].ix[data['E'] == 0.0] = 1e-256
        # put the evalues into a more reasonable range (thanks Richard!)
        data['E'] = -np.log10(data['E'])
        return preprocessing.scale(data)

    def fit(self, training_alignments):
        prepped = self.preprocess(training_alignments)
        self.classifier.fit(prepped)
        joblib.dump(self.classifier, self.model_fn)
        self.trained = True

    def predict(self, alignments):
        if not self.trained:
            self.classifier = joblib.load(self.model_fn)
        results = self.classifier.predict(alignments)
        return results == 1.0

    def fit_task(self):
        pass

    def predict_task(self):
        pass

