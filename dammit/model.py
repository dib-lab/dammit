#/usr/bin/env python
from __future__ import print_function

import numpy as np
import pandas as pd
from sklearn.externals import joblib
from sklearn import svm
from sklearn import preprocessing

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task

from . import common
from . import parsers
from . import tasks

class CRBL(object):

    def __init__(self, training_pep_fn, transcriptome_fn, database_pep_fn,
                 model_fn, classifier=svm.OneClassSVM, 
                 params={'nu':0.1, 'kernel':'rbf', 'gamma':0.1},
                 n_threads=1,
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
        self.train_rbh_fn = '{0}.rbhx.{1}.csv'.format(self.training_pep_fn,
                                                      self.database_pep_fn)
        self.crbl_fn = '{0}.crbl.{1}.csv'.format(self.transcriptome_fn,
                                                 self.database_pep_fn)

        self.model_fn = model_fn

        self.trained = False
        self.features = ['q_aln_len', 'E']
        self.classifier = classifier
        self.params = params
        self.model = classifier(**params)

        self.n_threads = n_threads
        self.lastal_cfg = lastal_cfg
        self.lastdb_cfg = lastdb_cfg

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

        CRBL.best_hits(aln_df_A)
        CRBL.best_hits(aln_df_B)

        # Join between subject A and query B
        rbh_df = pd.merge(aln_df_A, aln_df_B, how='inner', 
                          left_on='s_name', right_on='q_name')
        # Select those where query A is the same as subject B
        rbh_df = rbh_df[(rbh_df['q_name_x'] == rbh_df['s_name_y'])]

        # Drop extra columns and rename for clarity
        del rbh_df['s_name_x']
        del rbh_df['s_name_y']
        rbh_df.rename(columns={'q_name_x': 'q_name', 
                               'q_name_y': 's_name',
                               'q_aln_len_x': 'q_aln_len',
                               'E_x': 'E'},
                      inplace=True)
        
        return rbh_df

    def preprocess(self, alignments):

        # Subset out the features
        data = alignments.dropna()[self.features].astype(float)

        # e-values of 0.0 mess things up, set them to something really low
        data['E'].ix[data['E'] == 0.0] = 1e-256
        # put the evalues into a more reasonable range (thanks Richard!)
        data['E'] = -np.log10(data['E'])
        data = data.as_matrix()
        return preprocessing.scale(data)

    def fit(self, training_alignments):
        prepped = self.preprocess(training_alignments)
        print(type(prepped), type(prepped[0]), type(prepped[0,0]), prepped.ndim)
        for row in prepped:
            for item in row:
                if not isinstance(item, float):
                    print('non-float:', item, type(item))
        print(self.model)
        self.model.fit(prepped)
        print('fit!')
        joblib.dump(self.model, self.model_fn)
        self.trained = True

    def predict(self, alignments):
        if not self.trained:
            self.model = joblib.load(self.model_fn)
        prepped = self.preprocess(alignments)
        results = self.model.predict(prepped)
        return results == 1.0

    def format_training_pep_task(self):
        return tasks.get_lastdb_task(self.training_pep_fn,
                                     self.training_pep_fn,
                                     self.lastdb_cfg, 
                                     prot=True)

    def format_database_pep_task(self):
        return tasks.get_lastdb_task(self.database_pep_fn,
                                     self.database_pep_fn,
                                     self.lastdb_cfg,
                                     prot=True)

    def align_training_pep_task(self):
        return tasks.get_lastal_task(self.training_pep_fn,
                                     self.database_pep_fn,
                                     self.train_x_db_fn,
                                     False, self.n_threads,
                                     self.lastal_cfg)

    def align_database_pep_task(self):
        return tasks.get_lastal_task(self.database_pep_fn,
                                     self.training_pep_fn,
                                     self.db_x_train_fn,
                                     False, self.n_threads,
                                     self.lastal_cfg)

    def align_transcriptome_task(self):
        return tasks.get_lastal_task(self.transcriptome_fn,
                                     self.database_pep_fn,
                                     self.transcriptome_x_db_fn,
                                     False, self.n_threads,
                                     self.lastal_cfg)

    @tasks.create_task_object
    def training_reciprocals_task(self):

        def cmd():
            training_df = pd.concat([df for df in 
                                     parsers.maf_to_df_iter(self.train_x_db_fn)])
            db_df = pd.concat([df for df in
                               parsers.maf_to_df_iter(self.db_x_train_fn)])
            best_hits = CRBL.reciprocal_best_hits(training_df, db_df)
            best_hits.to_csv(self.train_rbh_fn)

        return {'name': '[CRBL]training_reciprocals:' + self.train_rbh_fn,
                'title': title_with_actions,
                'actions': [(cmd, [])],
                'file_dep': [self.train_x_db_fn, self.db_x_train_fn],
                'targets': [self.train_rbh_fn],
                'clean': [clean_targets]}

    @tasks.create_task_object
    def fit_task(self):
        
        def cmd():
            alignments = pd.read_csv(self.train_rbh_fn)
            self.fit(alignments)

        def clean():
            self.model = self.classifier(**self.params)
            self.trained = False

        return {'name': '[CRBL]fit:' + self.model_fn,
                'title': title_with_actions,
                'actions': [(cmd, [])],
                'file_dep': [self.train_rbh_fn],
                'targets': [self.model_fn],
                'clean': [clean_targets, (clean, [])]}

    @tasks.create_task_object
    def predict_task(self):
        
        def cmd():
            results = []
            for group in parsers.maf_to_df_iter(self.transcriptome_x_db_fn):
                mask = self.predict(group)
                results.append(group.ix[mask])
            pd.concat(results).to_csv(self.crbl_fn)

        return {'name': '[CRBL]predict:' + self.crbl_fn,
                'title': title_with_actions,
                'actions': [(cmd, [])],
                'file_dep': [self.model_fn, self.transcriptome_x_db_fn],
                'targets': [self.crbl_fn],
                'clean': [clean_targets]}

    def tasks(self):
        yield self.format_training_pep_task()
        yield self.format_database_pep_task()
        yield self.align_training_pep_task()
        yield self.align_database_pep_task()
        yield self.align_transcriptome_task()
        yield self.training_reciprocals_task()
        yield self.fit_task()
        yield self.predict_task()
