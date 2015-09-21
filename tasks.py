#!/usr/bin/env python
from __future__ import print_function

from itertools import izip
import json
import os
import pprint
import re
from shutil import rmtree
import shutil
import sys

from doit.tools import run_once, create_folder, title_with_actions, LongRunning
from doit.task import clean_targets, dict_to_task

from bioservices import UniProt
import jinja2
import pandas as pd
#import screed
from khmer import HLLCounter, ReadParser

def print_tasks(tasks):
    for task in tasks:
        print('-----\n', task)
        pprint.pprint(task.__dict__)

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
def get_uniprot_query_task(query, target_fn, fmt='fasta'):

    def func():
        u = UniProt()
        res = u.search(query, frmt=fmt)
        with open(target_fn, 'wb') as fp:
            fp.write(res)

    name = 'uniprot_query:' + query

    return {'name': name,
            'title': title_with_actions,
            'actions': [(func, [])],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [run_once] }

@create_task_object
def get_truncate_fasta_header_task(fasta_fn, length=500):
    
    def func():
        tmp_fn = fasta_fn + '.tmp'
        with open(tmp_fn, 'wb') as fp:
            for r in screed.open(fasta_fn):
                try:
                    name = r.name.split()[0]
                except IndexError:
                    name = r.name[:length]
                fp.write('>{}\n{}\n'.format(name, r.sequence))

        shutil.move(tmp_fn, fasta_fn)
    
    name = 'truncate_fasta_header:' + os.path.basename(fasta_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [(func, [])],
            'file_dep': [fasta_fn],
            'uptodate': [run_once]}

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
def get_blast_format_task(db_fn, db_out_fn, db_type):
    assert db_type in ['nucl', 'prot']

    cmd = 'makeblastdb -in {db_fn} -dbtype {db_type} -out {db_out_fn}'.format(**locals())

    target_fn = ''
    if db_type == 'nucl':
        target_fn = db_out_fn + '.nhr'
    else:
        target_fn = db_out_fn + '.phr'

    name = 'makeblastdb:' + os.path.basename(db_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [LongRunning(cmd), 'touch '+db_out_fn],
            'targets': [db_out_fn],
            'file_dep': [db_fn],
            'clean': [clean_targets, 'rm -f {target_fn}.*'.format(**locals())] }

def get_blast_task(query, db, prog, out_fn, n_threads, blast_dir, blast_cfg):
    assert prog in ['blastp', 'blastx', 'blastn', 'tblastn', 'tblastx']
    name = 'blast:' + os.path.basename(out_fn)

    params = blast_cfg['params']
    evalue = blast_cfg['evalue']
    exc = os.path.join(blast_dir, prog)

    cmd = '{exc} -query {query} -db {db} -num_threads {n_threads} '\
          '-evalue {evalue} -outfmt 6 {params} -o {out_fn}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [out_fn],
            'file_dep': [db],
            'clean': [clean_targets]}

@create_task_object
def link_file_task(src, dst=''):
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
def get_bowtie2_build_task(input_fn, db_basename, bowtie2_cfg):

    extra_args = bowtie2_cfg['extra_args']
    cmd = 'bowtie2-build {extra_args} {input_fn} {db_basename}'.format(**locals())
    
    targets = [db_basename+ext for ext in \
                ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']]
    targets.append(db_basename)

    name = 'bowtie2_build:' + os.path.basename(db_basename)
    return {'title': title_with_actions,
            'name': db_basename,
            'actions': [cmd, 'touch {db_basename}'.format(**locals())],
            'targets': targets,
            'file_dep': [input_fn],
            'clean': [clean_targets] }

@create_task_object
def get_bowtie2_align_task(db_basename, target_fn, bowtie2_cfg, n_threads, left_fn='', right_fn='', singleton_fn='',
                        read_fmt='-q', samtools_convert=True,
                        encoding='phred33'):

    assert read_fmt in ['-q', '-f']
    assert encoding in ['phred33', 'phred64']
    encoding = '--' + encoding
    if (left_fn or right_fn):
        assert (left_fn and right_fn)
    extra_args = bowtie2_cfg['extra_args']
    cmd = 'bowtie2 -p {n_threads} {extra_args} {encoding} {read_fmt} -x {db_basename} '.format(**locals())
    
    file_dep = [db_basename]
    targets = []

    name = 'bowtie2_align:' + ''.join('+' + os.path.basename(fn) if fn \
                                      else os.path.basename(fn) \
                                      for fn in [left_fn, right_fn, singleton_fn, db_basename])

    if left_fn:
        file_dep.extend([left_fn, right_fn])
        left_fn = '-1 ' + left_fn
        right_fn = '-2 ' + right_fn
    if singleton_fn:
        file_dep.append(singleton_fn)
        singleton_fn = '-U ' + singleton_fn
    if samtools_convert:
        targets.append(target_fn + '.bam')
        target_fn = ' | samtools view -Sb - > {target_fn}.bam'.format(**locals())
    else:
        targets.append(target_fn)
        target_fn = '-S ' + target_fn

    cmd = cmd + '{left_fn} {right_fn} {singleton_fn} {target_fn}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': targets,
            'file_dep': file_dep,
            'clean': [clean_targets] }

@create_task_object
def get_samtools_sort_task(bam_fn):
    
    cmd = 'samtools sort -n {bam_fn} {bam_fn}.sorted'.format(**locals())

    name = 'samtools_sort:' + os.path.basename(bam_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [bam_fn + '.sorted.bam'],
            'file_dep': [bam_fn],
            'clean': [clean_targets] }
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
def get_busco_task(input_filename, output_dir, busco_db_dir, input_type,
                   n_threads, busco_dir, busco_cfg):
    
    name = 'busco:' + os.path.basename(input_filename) + '-' + os.path.basename(busco_db_dir)

    assert input_type in ['genome', 'OGS', 'trans']
    exc = os.path.join(busco_dir, 'BUSCO_v1.1b1.py')

    cmd = 'python3 {exc} -in {in_fn} -o {out_dir} -l {db_dir} '\
            '-m {in_type} -c {n_threads}'.format(busco_path=busco_path, 
            in_fn=input_filename, out_dir=output_dir, db_dir=busco_db_dir, 
            in_type=input_type, n_threads=n_threads)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': ['run_' + output_dir, 
                        os.path.join('run_' + output_dir, 'short_summary_' + output_dir.rstrip('/'))],
            'file_dep': [input_filename],
            'uptodate': [run_once],
            'clean': [(clean_folder, ['run_' + output_dir])]}

@create_task_object
def get_cmpress_task(db_filename, infernal_dir, infernal_cfg):

    exc = os.path.join(infernal_dir, 'cmpress')
    cmd = '{exc} {db_filename}'.format(**locals())

    return {'name': 'cmpress:' + os.path.basename(db_filename),
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.i1f', '.i1i', '.i1m', '.i1p']],
            'file_dep': [db_filename],
            'clean': [clean_targets]}

@create_task_object
def get_cmscan_task(input_filename, output_filename, db_filename, 
                    n_threads, infernal_dir, infernal_cfg):
    
    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)
    
    exc = os.path.join(infernal_dir, 'cmscan')
    cmd = '{exc} --cpu {n_threads} --cut_ga --rfam --nohmmonly --tblout {output_filename}.tbl'\
          ' {db_filename} {input_filename} > {output_filename}.cmscan'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename + 'i3p'],
            'targets': [output_filename + '.tbl', output_filename + '.cmscan'],
            'clean': [clean_targets]}

@create_task_object
def get_hmmpress_task(db_filename, hmmer_dir, hmmer_cfg):
    
    name = 'hmmpress:' + os.path.basename(db_filename)
    exc = os.path.join(hmmer_dir, 'hmmpress')
    cmd = '{exc} {db_filename}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.h3f', '.h3i', '.h3m', '.h3p']],
            'file_dep': [db_filename],
            'clean': [clean_targets]}

@create_task_object
def get_hmmscan_task(input_filename, output_filename, db_filename, 
                     n_threads, hmmer_dir, hmmer_cfg):

    name = 'hmmscan:' + os.path.basename(input_filename) + '.x.' + \
                os.path.basename(db_filename)
    
    exc = os.path.join(hmmer_dir, 'hmmscan')
    cmd = '{exc} --cpu {n_threads} --domtblout {output_filename} \
          {db_filename} {input_filename}'.format(**locals())

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename+'.h3p'],
            'targets': [output_filename],
            'clean': [clean_targets]}

@create_task_object
def get_transdecoder_orf_task(input_filename, transdecoder_dir, transdecoder_cfg):

    name = 'TransDecoder.LongOrfs:' + os.path.basename(input_filename)

    min_prot_len = transdecoder_cfg['min_prot_len']
    exc = os.path.join(transdecoder_dir, 'TransDecoder.LongOrfs')
    cmd = '{exc} -t {input_filename} -m {min_prot_len}'.format(**locals())

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'targets': [input_filename + '.transdecoder_dir/longest_orfs.pep'],
            'clean': [(clean_folder, [input_filename + '.transdecoder_dir'])]}

# TransDecoder.Predict -t lamp10.fasta --retain_pfam_hits lamp10.fasta.pfam-A.out
@create_task_object
def get_transdecoder_predict_task(input_filename, db_filename, n_threads, 
                                  transdecoder_dir, transdecoder_cfg,):

    name = 'TransDecoder.Predict:' + os.path.basename(input_filename)

    orf_cutoff = transdecoder_cfg['orf_cutoff']
    exc = os.path.join(transdecoder_dir, 'TransDecoder.Predict')
    cmd = '{exc} -t {input_filename} --retain_pfam_hits {db_filename} \
            --retain_long_orfs {orf_cutoff} --cpu {n_threads}'.format(**locals())
    
    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, input_filename + '.transdecoder_dir/longest_orfs.pep', db_filename],
            'targets': [input_filename + ext for ext in ['.bed', '.cds', '.pep', '.gff3', '.mRNA']],
            'clean': [clean_targets, (clean_folder, [input_filename + '.transdecoder_dir'])]}

@create_task_object
def get_transcriptome_stats_task(transcriptome, out_dir):

    name = 'transcriptome_stats:' + os.path.basename(transcriptome)
    target = os.path.join(out_dir, transcriptome + '.stats')
    K = 25
    
    def parse(fn):
        hll = HLLCounter(.01, K)
        lens = []
        names = []
        for contig in ReadParser(fn):
            lens.append(len(contig.sequence))
            names.append(contig.name)
            hll.consume_string(contig.sequence)
        S = pd.Series(lens, index=names)
        S.sort()
        return S, hll.estimate_cardinality()

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
        lens, uniq_kmers = parse(transcriptome)
        
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
                 'redundancy': redundancy}
        
        with open(target, 'wb') as fp:
            json.dump(stats, fp, indent=4)

    return {'name': name,
            'title': title_with_actions,
            'actions': [(cmd, [])],
            'file_dep': [transcriptome],
            'targets': [target],
            'clean': [clean_targets]}
