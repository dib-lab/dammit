import os
from dammit.meta import __path__
from dammit.meta import get_databases

databases_d = get_databases()

rule download_and_gunzip:
    output: os.path.join(config["db_dir"], '{database}.{file_type}')
    params:
        url = lambda wildcards: databases_d[wildcards.database]['url'],
        md5 = lambda wildcards: databases_d[wildcards.database].get('md5', False),
        metalink = lambda wildcards: databases_d[wildcards.database].get('metalink', False),
        fileformat = lambda wildcards: databases_d[wildcards.database]['fileformat'],
        folder = lambda wildcards: databases_d[wildcards.database].get('folder', False)
    log: os.path.join(config["db_dir"], '{database}.{file_type}.log')
    wildcard_constraints:
        file_type = "hmm|cm|fasta|txt|done"
    script:
        f'file://{__path__}/wrappers/download/wrapper.py'

rule lastdb:
    input:
        os.path.join(config["db_dir"], "{database}.{file_type}")
    output:
        os.path.join(config["db_dir"], "{database}.{file_type}.prj"),
    params:
        protein_input =  lambda w: databases_d[w.database].get('db_type', False),
        extra = config["lastdb"]["params"].get('extra', " -w3 ")
    wildcard_constraints:
        file_type = "fasta|txt"
    log:
        "logs/{database}.{file_type}_lastdb.log"
    conda: f"file://{__path__}/wrappers/last/environment.yml"
    script: f"file://{__path__}/wrappers/last/lastdb.wrapper.py"

rule hmmpress:
    input:
        os.path.join(config["db_dir"], "{database}.hmm")
    output:
        os.path.join(config["db_dir"], "{database}.hmm.h3f"),
        os.path.join(config["db_dir"], "{database}.hmm.h3i"),
        os.path.join(config["db_dir"], "{database}.hmm.h3m"),
        os.path.join(config["db_dir"], "{database}.hmm.h3p")
    log:
        "logs/{database}_hmmpress.log"
    params:
        extra=config["hmmpress"]["params"].get("extra", ""),
    threads: 4
    conda: f"file://{__path__}/wrappers/hmmer/environment.yml"
    script: f"file://{__path__}/wrappers/hmmer/hmmpress.wrapper.py"

#rule hmmbuild:
#    input:
#        os.path.join(config["db_dir"], "{database}.sto")
#    output:
#        os.path.join(config["db_dir"], "{database}.hmm")
#    log:
#        "logs/{database}-hmmbuild.log"
#    params:
#        extra=config["hmmbuild"]["params"].get("extra", ""),
#    threads: 4
#    conda: f"file://{__path__}/wrappers/hmmer/environment.yml"
#    script: f"file://{__path__}/wrappers/hmmer/hmmbuild.wrapper.py"

rule infernal_cmpress:
    input:
        os.path.join(config["db_dir"], "{database}.cm")
    output:
        os.path.join(config["db_dir"],"{database}.cm.i1i"),
        os.path.join(config["db_dir"],"{database}.cm.i1f"),
        os.path.join(config["db_dir"],"{database}.cm.i1m"),
        os.path.join(config["db_dir"],"{database}.cm.i1p")
    log:
        "logs/cmpress_{database}.log"
    params:
        extra=config["cmpress"]["params"].get("extra", ""),
    conda: f"file://{__path__}/wrappers/infernal/environment.yml"
    script: f"file://{__path__}/wrappers/infernal/cmpress.wrapper.py"
