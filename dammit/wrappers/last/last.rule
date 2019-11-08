
import os

DATA_DIR="../../tests/test-data/"

## not working atm .. will figure it out later. these tests work in my snakemake-wrappers repo
# for now, test with: snakemake -s last.rule --use-conda

rule test_last:
    input: "test-tr-x-prot.maf", "test-transcript.maf"

rule lastdb_transcript:
    input:
        os.path.join(DATA_DIR, "test-transcript.fa")
    output:
        "test-transcript.fa.prj",
    params:
        protein_input = False,
        extra = ""
    log:
        "logs/lastdb/test-transcript.log"
    conda: "environment.yml"
    script: "lastdb.wrapper.py"

rule lastdb_protein:
    input:
        os.path.join(DATA_DIR, "test-protein.fa")
    output:
        "test-protein.fa.prj",
    params:
        protein_input = True,
        extra = ""
    log:
        "logs/lastdb/test-protein.log"
    conda: "environment.yml"
    script: "lastdb.wrapper.py"

rule lastal_nucl_x_nucl:
    input:
        data   = os.path.join(DATA_DIR, "test-transcript.fa"),
        lastdb ="test-transcript.fa.prj"
    output:
        maf = "test-transcript.maf"
    params:
        extra=""
    log:
        "logs/lastal/test.log"
    threads: 8
    conda: "environment.yml"
    script: "lastal.wrapper.py"

rule lastal_nucl_x_prot:
    input:
        data   = os.path.join(DATA_DIR, "test-transcript.fa"),
        lastdb ="test-protein.fa.prj"
    output:
        maf = "test-tr-x-prot.maf"
    params:
        extra=""
    log:
        "logs/lastal/test.log"
    threads: 8
    conda: "environment.yml"
    script: "lastal.wrapper.py"
