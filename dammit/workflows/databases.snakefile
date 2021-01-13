import os
from dammit.meta import __path__, __wrappers__
from dammit.config import DATABASES

THREADS_PER_TASK = config['max_threads_per_task']


rule download_and_gunzip:
    message:
        """
        Download and decompress a remote file.
        """
    output: os.path.join(config["database_dir"], '{database}.{file_type}')
    params:
        url = lambda wildcards: DATABASES[wildcards.database]['url'],
        md5 = lambda wildcards: DATABASES[wildcards.database].get('md5', False),
        metalink = lambda wildcards: DATABASES[wildcards.database].get('metalink', False),
        fileformat = lambda wildcards: DATABASES[wildcards.database].get('fileformat', False),
        folder = lambda wildcards: DATABASES[wildcards.database].get('folder', False)
    log: os.path.join(logs_dir, '{database}.{file_type}.log')
    wildcard_constraints:
        file_type = "hmm|cm|fasta|txt|ini|done"
    threads: 1
    script: f'file://{__wrappers__}/download/wrapper.py'


rule lastdb:
    message:
        """
        Prepare a protein FASTA for use as a database by the LAST aligner.
        """
    input:
        os.path.join(config["database_dir"], "{database}.{file_type}")
    output:
        os.path.join(config["database_dir"], "{database}.{file_type}.prj"),
    params:
        protein_input =  lambda w: DATABASES[w.database].get('db_type', False),
        extra = config["lastdb"]["params"].get('extra', " -w3 ")
    wildcard_constraints:
        file_type = "fasta|txt"
    log:
        os.path.join(results_dir, "{database}.{file_type}_lastdb.log")
    threads: THREADS_PER_TASK
    wrapper: f"file://{__wrappers__}/last/lastdb.wrapper.py"


rule hmmpress:
    message:
        """
        Prepare a collection of profile hidden markov models (Pfam-A)
        for use by HMMER.
        """
    input:
        os.path.join(config["database_dir"], "{database}.hmm")
    output:
        os.path.join(config["database_dir"], "{database}.hmm.h3f"),
        os.path.join(config["database_dir"], "{database}.hmm.h3i"),
        os.path.join(config["database_dir"], "{database}.hmm.h3m"),
        os.path.join(config["database_dir"], "{database}.hmm.h3p")
    log:
        os.path.join(results_dir, "{database}_hmmpress.log")
    params:
        extra=config["hmmpress"]["params"].get("extra", ""),
    threads: 1
    wrapper: f'file://{__wrappers__}/hmmer/hmmpress.wrapper.py'


rule infernal_cmpress:
    message:
        """
        Prepare a collection of covariance models (Rfam) for use by Infernal.
        """
    input:
        os.path.join(config["database_dir"], "{database}.cm")
    output:
        os.path.join(config["database_dir"],"{database}.cm.i1i"),
        os.path.join(config["database_dir"],"{database}.cm.i1f"),
        os.path.join(config["database_dir"],"{database}.cm.i1m"),
        os.path.join(config["database_dir"],"{database}.cm.i1p")
    log:
        os.path.join(results_dir, "cmpress_{database}.log")
    params:
        extra=config["cmpress"]["params"].get("extra", ""),
    threads: 1
    wrapper: f'file://{__wrappers__}/infernal/cmpress.wrapper.py'
