import os
from dammit.meta import __path__

configfile: '../../databases.json'
DATABASES_DIR= f"{__path__}/tests/databases/"

rule download_databases:
    input: 'Pfam-A.hmm', 'Rfam.cm', 'OrthoDB.fasta', 'uniprot_sprot.fasta'


rule download_and_gunzip:
    output: '{database}.{file_type}'
    params:
        url = lambda wildcards: config[wildcards.database]['url'],
        md5 = lambda wildcards: config[wildcards.database].get('md5', False),
        metalink = lambda wildcards: config[wildcards.database].get('metalink', False)
    wrapper:
        f'file://{__path__}/wrappers/download'

rule lastdb:
    input:
        os.path.join(DATABASES_DIR, "{database}.fa")
    output:
        os.path.join(DATABASES_DIR, "{database}.fa.prj"),
    params:
        protein_input =  lambda w: config[w.database].get('db_type', False),
        extra = config["lastdb"]["params"].get('extra', " -w3 ")
    log:
        "logs/{database}_lastdb.log"
    conda: f"file://{__path__}/wrappers/last/environment.yml"
    wrapper: f"file://{__path__}/wrappers/last/lastdb.wrapper.py"

rule hmmpress:
    input:
        os.path.join(DATABASES_DIR, "{database}.hmm")
    output:
        os.path.join(DATABASES_DIR, "{database}.hmm.h3f"),
        os.path.join(DATABASES_DIR, "{database}.hmm.h3i"),
        os.path.join(DATABASES_DIR, "{database}.hmm.h3m"),
        os.path.join(DATABASES_DIR, "{database}.hmm.h3p")
    log:
        "logs/{database}_hmmpress.log"
    params:
        extra=config["hmmpress"]["params"].get("extra", ""),
    threads: 4
    conda: f"file://{__path__}/wrappers/hmmer/environment.yml"
    wrapper: f"file://{__path__}/wrappers/hmmer/hmmpress.wrapper.py"

rule hmmbuild:
    input:
        os.path.join(DATABASES_DIR, "{database}.sto")
    output:
        os.path.join(DATABASES_DIR, "{database}.hmm")
    log:
        "logs/{database}-hmmbuild.log"
    params:
        extra=config["hmmbuild"]["params"].get("extra", ""),
    threads: 4
    conda:
        f"file://{__path__}/wrappers/hmmer/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmbuild.wrapper.py"

rule infernal_cmpress:
    input:
        os.path.join(DATA_DIR, "transcriptome.cm")
    output:
        os.path.join(DATA_DIR,"{transcriptome}.cm.i1i"),
        os.path.join(DATA_DIR,"{transcriptome}.cm.i1f"),
        os.path.join(DATA_DIR,"{transcriptome}.cm.i1m"),
        os.path.join(DATA_DIR,"{transcriptome}.cm.i1p")
    log:
        "logs/cmpress.log"
    params:
        extra=config["cmpress"]["params"].get("extra", ""),
    conda: f"file://{__path__}/wrappers/infernal/environment.yml"
    wrapper: f"file://{__path__}/wrappers/infernal/cmpress.wrapper.py"
