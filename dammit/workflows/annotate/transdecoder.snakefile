import os
from dammit.meta import __path__

DATA_DIR = f"{__path__}/tests/test-data"

# snakemake -s transdecoder.snakefile --use-conda 
rule test_transdecoder:
    input: "test.fa.transdecoder.bed"
    
rule transdecoder_longorfs:
    input:
        fasta = os.path.join(DATA_DIR, "test.fa.gz")
    output:
        "{assembly}.transdecoder_dir/longest_orfs.pep"
    log:
        'logs/transdecoder/{assembly}-longorfs.log'
    params:
        extra=""
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-longorfs.wrapper.py"

#rule hmmer_pfam:
#    input:
#    output:
#    log:
#    params:
#    threads:
#    conda: f'file://{__path__}/wrappers/hmmer/environment.yml'
#    wrapper:
#        f'file://{__path__}/wrappers/hmmer'

rule transdecoder_predict:
    input:
        fasta = os.path.join(DATA_DIR, "{assembly}.gz"),
        longorfs = "{assembly}.transdecoder_dir/longest_orfs.pep"
    output:
        "{assembly}.transdecoder.bed",
        "{assembly}.transdecoder.cds",
        "{assembly}.transdecoder.pep",
        "{assembly}.transdecoder.gff3"
    log:
        'logs/transdecoder/{assembly}-predict.log'
    params:
        extra=""
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-predict.wrapper.py"
