#!/usr/bin/env python
from __future__ import print_function

from itertools import count, izip
import json
import logging
import os
import pprint
import re
from shutil import rmtree
import shutil
import sys

from doit.tools import run_once, create_folder, title_with_actions, LongRunning
from doit.task import clean_targets, dict_to_task

import pandas as pd
from khmer import HLLCounter, ReadParser

from common import which
from . import parsers
from .hits import BestHits
#from . import taxonomy
import gff


def task_str(task):
    return pprint.pformat(task.__dict__)


def print_tasks(tasks, logger=None, level=logging.DEBUG):
    for task in tasks:
        if logger:
            logger.log(level, task_str(task))
        else:
            print(task, file=sys.stderr)


def clean_folder(target):
    try:
        rmtree(target)
    except OSError:
        pass


seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]


def create_task_object(task_dict_func):
    '''Wrapper to decorate functions returning pydoit
    Task dictionaries and have them return pydoit Task
    objects
    '''
    def d_to_t(*args, **kwargs):
        ret_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(ret_dict)
    return d_to_t


@create_task_object
def get_group_task(group_name, tasks):

    return {'name': group_name,
            'actions': None,
            'task_dep': [t.name for t in tasks]}


@create_task_object
def get_download_task(url, target_fn, label='default'):

    cmd = 'curl -o {target_fn} {url}'.format(**locals())
    name = '_'.join(['download_gunzip', target_fn, label])

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}


@create_task_object
def get_download_and_gunzip_task(url, target_fn):
    cmd = 'curl {url} | gunzip -c > {target_fn}'.format(**locals())

    name = 'download_and_gunzip:' + os.path.basename(target_fn)

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}


@create_task_object
def get_download_and_untar_task(url, target_dir, label):

    cmd1 = 'mkdir -p {target_dir}; curl {url} | tar -xz -C {target_dir}'.format(**locals())
    name = 'download_and_untar:'  + \
            os.path.basename(os.path.dirname(target_dir)) + '-' + label
    done = os.path.join(target_dir, name + '.done')
    cmd2 = 'touch {done}'.format(done=done)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd1, cmd2],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [True]}


@create_task_object
def get_create_folder_task(folder):

    name = 'create_folder:{folder}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [(create_folder, [folder])],
            'targets': [folder],
            'uptodate': [run_once],
            'clean': [clean_targets] }


@create_task_object
def get_sanitize_fasta_task(input_fn, output_fn):

    import re
    name = 'sanitize_fasta:{0}'.format(os.path.basename(input_fn))

    def fix():
        bad = r'["`/|=]+'
        with open(output_fn, 'wb') as fp:
            for record in ReadParser(input_fn):
                header = re.sub(bad, r'.', record.name)
                fp.write('>{0}\n{1}\n'.format(header, record.sequence))

    return {'name': name,
            'title': title_with_actions,
            'actions': [fix],
            'file_dep': [input_fn],
            'targets': [output_fn],
            'clean': [clean_targets]}


@create_task_object
def get_rename_transcriptome_task(transcriptome_fn, output_fn, names_fn, 
                                  transcript_basename, split_regex=None):

    import re
    name = os.path.basename(transcriptome_fn)

    if split_regex is None:
        counter = count()
        header_func = lambda name: '{0}_{1}'.format(transcript_basename, next(counter))
    else:
        def header_func(header):
            results = re.search(split_regex, header).groupdict()
            try:
                header = results['name']
            except KeyError as err:
                err.message = 'Header regex should have a name field!'
                raise
            return header

    def fix():
        names = []
        with open(output_fn, 'wb') as fp:
            for record in ReadParser(transcriptome_fn):
                header = header_func(record.name)
                fp.write('>{0}\n{1}\n'.format(header, record.sequence))
                names.append((record.name, header))
        pd.DataFrame(names, columns=['original', 'renamed']).to_csv(names_fn,
                                                                    index=False)

    return {'name': name,
            'title': title_with_actions,
            'actions': [fix],
            'targets': [output_fn, names_fn],
            'file_dep': [transcriptome_fn],
            'clean': [clean_targets]} 


@create_task_object
def get_transeq_task(input_fn, output_fn, clean=True, frame=6):

    name = 'transeq:{0}'.format(os.path.basename(input_fn))

    params = '-frame {frame}'.format(frame=frame)
    if clean:
        params = '{params} -clean'.format(params=params)

    cmd = 'transeq -sequence {input_fn} {params}'\
          ' -clean -outseq {output_fn}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}


@create_task_object
def get_crb_blast_task(query, target, output, cutoff, crb_blast_cfg,
                       n_threads):

    name = 'crb-blast:{0}.x.{1}'.format(query, target)
    exc = which('crb-blast')
    cmd = '{exc} --query {query} --target {target} --output {output} '\
          '--evalue {cutoff} --threads {n_threads}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output],
            'file_dep': [query, target],
            'clean': [clean_targets]}


@create_task_object
def get_lastdb_task(db_fn, db_out_prefix, lastdb_cfg, prot=True):
    '''Create a pydoit task to run lastdb.

    WARNING: This does not define a file_dep, to make sure it doesn't
    get executed when the dependency and targets already exist. This means
    that if a task acquires the database, it MUST be defined before the
    lastdb task.

    Args:
        db_fn (str): The FASTA file to format.
        db_out_prefix (str): Prefix for the database files.
        lastdb_cfg (dict): Config for the command. Shoud contain an entry
        named "params" storing a str.
        prot (bool): True if a protein FASTA, False otherwise.
    Returns:
        dict: A pydoit task.
    '''
    
    exc = which('lastdb')
    params = lastdb_cfg['params']
    if prot:
        params += ' -p'

    cmd = '{exc} {params} {db_out_prefix} {db_fn}'.format(**locals())

    name = 'lastdb:' + os.path.basename(db_out_prefix)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_out_prefix + ext \
                        for ext in \
                        ['.bck', '.des', '.prj', '.sds', '.ssp', '.suf', '.tis']],
            'uptodate': [True],
            'clean': [clean_targets]}


@create_task_object
def get_lastal_task(query, db, out_fn, translate, cutoff, n_threads, lastal_cfg):
    '''Create a pydoit task to run lastal

    Args:
        query (str): The file with the query sequences.
        db (str): The database file prefix.
        out_fn (str): Destination file for alignments.
        translate (bool): True if query is a nucleotide FASTA.
        n_threads (int): Number of threads to run with.
        lastal_cfg (dict): Config, must contain key params holding str.
    Returns:
        dict: A pydoit task.
    '''

    exc = which('lastal')
    params = lastal_cfg['params']

    if translate:
        params += ' -F' + str(lastal_cfg['frameshift'])
    if cutoff is not None:
        cutoff = 1.0 / cutoff
        params += ' -D' + str(cutoff)
    cmd = '{exc} {params} {db} {query} > {out_fn}'.format(**locals())

    name = 'lastal:' + os.path.join(out_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [out_fn],
            'file_dep': [db + '.sds'],
            'clean': [clean_targets]}


@create_task_object
def get_maf_best_hits_task(maf_fn, output_fn):

    hits_mgr = BestHits()

    def cmd():
        df = pd.concat([group for group in
                        parsers.maf_to_df_iter(maf_fn)])
        df = hits_mgr.best_hits(df)
        df.to_csv(output_fn, index=False)

    name = 'maf_best_hits:{0}-{1}'.format(maf_fn, output_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [maf_fn],
            'clean': [clean_targets]}



@create_task_object
def get_link_file_task(src, dst=''):
    ''' Soft-link file to the current directory
    '''
    cmd = 'ln -fs {src} {dst}'.format(src=src, dst=dst)
    return {'title': title_with_actions,
            'name': 'ln:' + os.path.basename(src) + ('-' + dst if dst else ''),
            'actions': [cmd],
            'file_dep': [src],
            'targets': [os.path.basename(src) if not dst else dst],
            'uptodate': [run_once],
            'clean': [clean_targets]}


@create_task_object
def get_cat_task(file_list, target_fn):

    cmd = 'cat {files} > {t}'.format(files=' '.join(file_list), t=target_fn)

    return {'title': title_with_actions,
            'name': 'cat:' + os.path.basename(target_fn),
            'actions': [cmd],
            'file_dep': file_list,
            'targets': [target_fn],
            'clean': [clean_targets]}


@create_task_object
def get_busco_task(input_filename, output_name, busco_db_dir, input_type,
                   n_threads, busco_cfg):
    
    name = 'busco:' + os.path.basename(input_filename) + '-' + os.path.basename(busco_db_dir)

    assert input_type in ['genome', 'OGS', 'trans']
    exc = which('BUSCO_v1.1b1.py')
    # BUSCO chokes on file paths as output names
    output_name = os.path.basename(output_name)

    cmd = 'python3 {exc} -in {input_filename} -f -o {output_name} -l {busco_db_dir} '\
            '-m {input_type} -c {n_threads}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'uptodate': [run_once],
            'clean': [(clean_folder, ['run_' + output_name])]}


@create_task_object
def get_cmpress_task(db_filename, infernal_cfg):

    exc = which('cmpress')
    cmd = '{exc} {db_filename}'.format(**locals())

    return {'name': 'cmpress:' + os.path.basename(db_filename),
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.i1f', '.i1i', '.i1m', '.i1p']],
            'uptodate': [True],
            'clean': [clean_targets]}


@create_task_object
def get_cmscan_task(input_filename, output_filename, db_filename, 
                    cutoff, n_threads, infernal_cfg):
    
    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)
    
    exc = which('cmscan')
    cmd = '{exc} --cpu {n_threads} --rfam --nohmmonly -E {cutoff}'\
          ' --tblout {output_filename} {db_filename} {input_filename}'\
          ' > {output_filename}.cmscan'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename + '.i1p'],
            'targets': [output_filename, output_filename + '.cmscan'],
            'clean': [clean_targets]}


@create_task_object
def get_hmmpress_task(db_filename, hmmer_cfg):
    
    name = 'hmmpress:' + os.path.basename(db_filename)
    exc = which('hmmpress')
    cmd = '{exc} {db_filename}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.h3f', '.h3i', '.h3m', '.h3p']],
            'uptodate': [True],
            'clean': [clean_targets]}


@create_task_object
def get_hmmscan_task(input_filename, output_filename, db_filename, 
                     cutoff, n_threads, hmmer_cfg):

    name = 'hmmscan:' + os.path.basename(input_filename) + '.x.' + \
                os.path.basename(db_filename)
    
    exc = which('hmmscan')
    stat = output_filename + '.out'
    cmd = '{exc} --cpu {n_threads} --domtblout {output_filename} -E {cutoff}'\
          ' -o {stat} {db_filename} {input_filename}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename+'.h3p'],
            'targets': [output_filename, stat],
            'clean': [clean_targets]}


@create_task_object
def get_transdecoder_orf_task(input_filename, transdecoder_cfg):

    name = 'TransDecoder.LongOrfs:' + os.path.basename(input_filename)

    min_prot_len = transdecoder_cfg['min_prot_len']
    exc = which('TransDecoder.LongOrfs')
    cmd = '{exc} -t {input_filename} -m {min_prot_len}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'targets': [input_filename + '.transdecoder_dir/longest_orfs.pep'],
            'clean': [(clean_folder, [input_filename + '.transdecoder_dir'])]}


@create_task_object
def get_remap_hmmer_task(hmmer_filename, remap_gff_filename, output_filename):

    name = 'remap_hmmer:{0}'.format(os.path.basename(hmmer_filename))

    def cmd():
        gff_df = pd.concat([group for group in
                            parsers.parse_gff3(remap_gff_filename)])
        hmmer_df = pd.concat([group for group in
                              parsers.hmmscan_to_df_iter(hmmer_filename)])

        merged_df = pd.merge(hmmer_df, gff_df, left_on='full_query_name', right_on='ID')

        hmmer_df['env_coord_from'] = (merged_df.start + \
                                      (3 * merged_df.env_coord_from)).astype(int)
        hmmer_df['env_coord_to'] = (merged_df.start + \
                                    (3 * merged_df.env_coord_to)).astype(int)
        assert(all(hmmer_df.env_coord_to <= merged_df.end))

        hmmer_df.to_csv(output_filename, header=True, index=False)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [hmmer_filename, remap_gff_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


# TransDecoder.Predict -t lamp10.fasta --retain_pfam_hits lamp10.fasta.pfam-A.out
@create_task_object
def get_transdecoder_predict_task(input_filename, db_filename, transdecoder_cfg):

    name = 'TransDecoder.Predict:' + os.path.basename(input_filename)

    orf_cutoff = transdecoder_cfg['orf_cutoff']
    exc = which('TransDecoder.Predict')
    cmd = '{exc} -t {input_filename} --retain_pfam_hits {db_filename} \
            --retain_long_orfs {orf_cutoff}'.format(**locals())
    
    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, 
                         input_filename + '.transdecoder_dir/longest_orfs.pep',
                         db_filename],
            'targets': [input_filename + '.transdecoder' + ext \
                        for ext in ['.bed', '.cds', '.pep', '.gff3', '.mRNA']],
            'clean': [clean_targets, 
                     (clean_folder, [input_filename + '.transdecoder_dir'])]}


@create_task_object
def get_maf_gff3_task(input_filename, output_filename, database):

    name = 'maf-gff3:' + os.path.basename(output_filename)

    def cmd():
        if input_filename.endswith('.csv') or input_filename.endswith('.tsv'):
            it = pd.read_csv(input_filename, chunksize=10000)
        else:
            it = parsers.maf_to_df_iter(input_filename)

        with open(output_filename, 'a') as fp:
            for group in it:
                gff_group = gff.maf_to_gff3_df(group, 'dammit.last', 
                                               database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@create_task_object
def get_crb_gff3_task(input_filename, output_filename, database):

    name = 'crbb-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in parsers.crb_to_df_iter(input_filename,
                                                remap=True):
                gff_group = gff.crb_to_gff3_df(group, 'dammit.crbb', database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@create_task_object
def get_hmmscan_gff3_task(input_filename, output_filename, database):

    name = 'hmmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in pd.read_csv(input_filename, chunksize=10000):
                print(group.head())
                gff_group = gff.hmmscan_to_gff3_df(group, 'dammit.hmmscan',
                                                   database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@create_task_object
def get_cmscan_gff3_task(input_filename, output_filename, database):

    name = 'cmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in parsers.cmscan_to_df_iter(input_filename):
                gff_group = gff.cmscan_to_gff3_df(group, 'dammit.cmscan',
                                                  database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@create_task_object
def get_gff3_merge_task(gff3_filenames, output_filename):

    name = 'gff3-merge:{0}'.format(os.path.basename(output_filename))

    merge_cmd = 'echo "{v}" > {out}; cat {f} | sed \'/^#/ d\''\
                ' | sort | sed \'/^$/d\' >> {out}'.format(v=gff.version_line(),
                                          f=' '.join(gff3_filenames),
                                          out=output_filename)
    return {'name': name,
            'title': title_with_actions,
            'actions': [merge_cmd],
            'file_dep': gff3_filenames,
            'targets': [output_filename],
            'clean': [clean_targets]}


@create_task_object
def get_annotate_fasta_task(transcriptome_fn, gff3_fn, output_fn):
    
    name = 'fasta-annotate:{0}'.format(output_fn)

    def annotate_fasta():
        annotations = pd.concat([g for g in \
                                 parsers.parse_gff3(gff3_fn)])
        with open(output_fn, 'wb') as fp:
            for n, record in enumerate(ReadParser(transcriptome_fn)):
                df = annotations.query('seqid == "{0}"'.format(record.name))
                annots = ['len={0}'.format(len(record.sequence))]
                #for seqid, sgroup in df.groupby('seqid'):
                #    for feature_type, fgroup in sgroup.groupby('feature_type'):
                for feature_type, fgroup in df.groupby('feature_type'):

                    if feature_type in ['translated_nucleotide_match',
                                        'protein_hmm_match',
                                        'RNA_sequence_secondary_structure']:

                        collapsed = ','.join(['{}:{}-{}'.format(row.Name.split(':dammit')[0], 
                                                                 int(row.start), 
                                                                 int(row.end)) \
                                        for _, row in fgroup.iterrows()])
                        if feature_type == 'translated_nucleotide_match':
                            key = 'homologies'
                        elif feature_type == 'protein_hmm_match':
                            key = 'hmm_matches'
                        else:
                            key = 'RNA_matches'
                        annots.append('{0}={1}'.format(key, collapsed))

                    elif feature_type in ['exon', 'CDS', 'gene',
                                          'five_prime_UTR', 'three_prime_UTR', 
                                          'mRNA']:
                        
                        collapsed = ','.join(['{}-{}'.format(int(row.start),
                                                             int(row.end)) \
                                        for _, row in fgroup.iterrows()])
                        annots.append('{0}={1}'.format(feature_type, collapsed))

                desc = '{0} {1}'.format(record.name, ' '.join(annots))

                fp.write('>{0}\n{1}\n'.format(desc.strip(), record.sequence))

    return {'name': name,
                'title': title_with_actions,
                'actions': [annotate_fasta],
                'file_dep': [transcriptome_fn, gff3_fn],
                'targets': [output_fn],
                'clean': [clean_targets]}


@create_task_object
def get_transcriptome_stats_task(transcriptome, output_fn):

    name = 'transcriptome_stats:' + os.path.basename(transcriptome)
    K = 25
    
    def parse(fn):
        hll = HLLCounter(.01, K)
        lens = []
        names = []
        gc_len = 0
        for contig in ReadParser(fn):
            lens.append(len(contig.sequence))
            names.append(contig.name)
            hll.consume_string(contig.sequence)
            gc_len += contig.sequence.count('C')
            gc_len += contig.sequence.count('G')
        S = pd.Series(lens, index=names)
        try:
            S.sort_values()
        except AttributeError:
            S.sort()
        gc_perc = float(gc_len) / S.sum()
        return S, hll.estimate_cardinality(), gc_perc

    def calc_NX(lens, X):
        N = lens.sum()
        threshold = (float(X) / 100.0) * N

        NXlen = 0
        NXpos = 0
        cms = 0
        for n, l in enumerate(lens):
            cms += l
            if cms >= threshold:
                NXlen = l
                NXpos = n
                break
        return NXlen, NXpos

    def cmd():
        lens, uniq_kmers, gc_perc = parse(transcriptome)
        
        exp_kmers = (lens - (K+1)).sum()
        redundancy = float(exp_kmers - uniq_kmers) / exp_kmers
        if redundancy < 0:
            redundancy = 0.0

        N50len, N50pos = calc_NX(lens, 50)
        stats = {'N': len(lens),
                 'sum': lens.sum(),
                 'min': lens.min(),
                 'max': lens.max(),
                 'med': lens.median(),
                 'mean': lens.mean(),
                 'N50len': N50len,
                 'N50pos': N50pos,
                 '25_mers': exp_kmers,
                 '25_mers_unique': uniq_kmers,
                 'redundancy': redundancy,
                 'GCperc': gc_perc}
        
        with open(output_fn, 'wb') as fp:
            json.dump(stats, fp, indent=4)

    return {'name': name,
            'title': title_with_actions,
            'actions': [(cmd, [])],
            'file_dep': [transcriptome],
            'targets': [output_fn],
            'clean': [clean_targets]}
