import os
from dammit.meta import __path__

DATA_DIR= f"{__path__}/tests/test-data/"

rule test_hmmer:
    input: "hmmscan/test-prot-tbl.txt", "hmmsearch/test-prot-tbl.txt"

rule hmmbuild_profile:
    input:
        os.path.join(DATA_DIR, "test-profile.sto")
    output:
        "test-profile.hmm"
    log:
        "logs/test-profile-hmmbuild.log"
    params:
        extra="",
    threads: 4
    conda:
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmbuild.wrapper.py"

rule hmmpress_profile:
    input:
        "test-profile.hmm"
    output:
        "test-profile.hmm.h3f",
        "test-profile.hmm.h3i",
        "test-profile.hmm.h3m",
        "test-profile.hmm.h3p"
    log:
        "logs/hmmpress.log"
    params:
        extra="",
    threads: 4
    conda:
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmpress.wrapper.py"

rule hmmscan_profile:
    input:
        fasta=os.path.join(DATA_DIR,"test-protein.fa"),
        profile="test-profile.hmm.h3f",
    output:
        # only one of these is required
        tblout="hmmscan/test-prot-tbl.txt", # save parseable table of per-sequence hits to file <f>
        domtblout="hmmscan/test-prot-domtbl.txt", # save parseable table of per-domain hits to file <f>
        pfamtblout="hmmscan/test-prot-pfamtbl.txt", # save table of hits and domains to file, in Pfam format <f>
        outfile="hmmscan/test-prot-out.txt", # Direct the main human-readable output to a file <f> instead of the default stdout.
    log:
        "logs/hmmscan.log"
    params:
        evalue_threshold=0.00001,
        # if bitscore threshold provided, hmmscan will use that instead
        #score_threshold=50,
        extra="",
    threads: 4
    conda:
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmscan.wrapper.py"

rule hmmsearch_profile:
    input:
        fasta=os.path.join(DATA_DIR,"test-protein.fa"),
        profile="test-profile.hmm.h3f",
    output:
        # only one of these is required
        tblout="hmmsearch/test-prot-tbl.txt", # save parseable table of per-sequence hits to file <f>
        domtblout="hmmsearch/test-prot-domtbl.txt", # save parseable table of per-domain hits to file <f>
        alignment_hits="hmmsearch/test-prot-alignment-hits.txt", # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
        outfile="hmmsearch/test-prot-out.txt", # Direct the main human-readable output to a file <f> instead of the default stdout. 
    log:
        "logs/hmmsearch.log"
    params:
        evalue_threshold=0.00001,
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra="",
    threads: 4
    conda:
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmsearch.wrapper.py"
