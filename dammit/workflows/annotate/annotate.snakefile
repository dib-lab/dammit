import os
from dammit.meta import __path__

DATA_DIR= f"{__path__}/tests/test-data/"
DATABASES_DIR= f"{__path__}/tests/databases/"
OUTDIR= f"{__path__}/tests/dammit-out/"

rule transdecoder_longorfs:
    input:
        fasta = os.path.join(DATA_DIR, "{transcriptome}.fa")
    output:
        "{transcriptome}.transdecoder_dir/longest_orfs.pep"
    log:
        'logs/transdecoder/{transcriptome}-longorfs.log'
    params:
        extra=config["transdecoder_longorfs"]["params"].get("extra", "-m 80 ")
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-longorfs.wrapper.py"

rule transdecoder_predict:
    input:
        fasta = os.path.join(DATA_DIR, "{transcriptome}.fa"),
        longorfs = os.path.join(DATABASES_DIR, "{transcriptome}.transdecoder_dir/longest_orfs.pep")
    output:
        os.path.join(DATABASES_DIR, "{transcriptome}.transdecoder.bed"),
        os.path.join(DATABASES_DIR, "{transcriptome}.transdecoder.cds"),
        os.path.join(DATABASES_DIR, "{transcriptome}.transdecoder.pep"),
        os.path.join(DATABASES_DIR, "{transcriptome}.transdecoder.gff3")
    log:
        'logs/transdecoder/{transcriptome}-predict.log'
    params:
        extra= config["transdecoder_predict"]["params"].get("extra", "")
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-predict.wrapper.py"

rule lastal:
    input:
        data   = os.path.join(DATA_DIR, "{transcriptome}.fa"),
        lastdb =os.path.join(DATABASES_DIR, "{transcriptome}.fa.prj")
    output:
        maf = os.path.join(OUTDIR, "{transcriptome}.maf")
    params:
        frameshift_cost = config["lastal"]["params"].get("frameshift_cost", 15),
        extra=config["lastal"]["params"].get("extra", ""),
    log:
        "logs/lastal/{transcriptome}.log"
    threads: 8
    conda: f"file://{__path__}/wrappers/last/environment.yml"
    script: f"file://{__path__}/wrappers/last/lastal.wrapper.py"

# probably want to switch to hmmsearch instead of hmmscan
rule hmmscan:
    input:
        fasta=os.path.join(DATA_DIR, "{transcriptome}.fa"),
        profile=os.path.join(DATABASES_DIR, "{database}.hmm.h3f"),
    output:
        # only one of these is required
        tblout=os.path.join(OUTDIR, "hmmscan/{transcriptome}-tbl.txt"), # save parseable table of per-sequence hits to file <f>
        domtblout=os.path.join(OUTDIR,"hmmscan/{transcriptome}-domtbl.txt"), # save parseable table of per-domain hits to file <f>
        pfamtblout=os.path.join(OUTDIR,"hmmscan/{transcriptome}-pfamtbl.txt"), # save table of hits and domains to file, in Pfam format <f>
        outfile=os.path.join(OUTDIR,"hmmscan/{transcriptome}-out.txt"), # Direct the main human-readable output to a file <f> instead of the default stdout.
    log:
        "logs/{transcriptome}_hmmscan.log"
    params:
        evalue_threshold=0.00001,
        # if bitscore threshold provided, hmmscan will use that instead
        #score_threshold=50,
        extra=config["hmmscan"]["params"].get("extra", ""),
    threads: 4
    conda:
        f"file://{__path__}/wrappers/hmmer/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/hmmer/hmmscan.wrapper.py"

rule hmmsearch:
    input:
        fasta=os.path.join(DATA_DIR, "{transcriptome}.fa"),
        profile=os.path.join(DATABASES_DIR, "{database}.hmm.h3f"),
    output:
        # only one of these is required
        tblout=os.path.join(OUTDIR, "hmmsearch/{transcriptome}-tbl.txt"), # save parseable table of per-sequence hits to file <f>
        domtblout=os.path.join(OUTDIR,"hmmsearch/{transcriptome}-domtbl.txt", # save parseable table of per-domain hits to file <f>
        alignment_hits=os.path.join(OUTDIR,"hmmsearch/{transcriptome}-alignment-hits.txt", # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
        outfile=os.path.join(OUTDIR,"hmmsearch/{transcriptome}-out.txt", # Direct the main human-readable output to a file <f> instead of the default stdout. 
    log:
        "logs/{transcriptome}_hmmsearch.log"
    params:
        evalue_threshold=0.00001,
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra=config["hmmsearch"]["params"].get("extra", ""),
    threads: 4
    conda: f"file://{__path__}/wrappers/hmmer/environment.yml"
    wrapper: f"file://{__path__}/wrappers/hmmer/hmmsearch.wrapper.py"

rule cmscan:
    input:
        fasta=  os.path.join(DATA_DIR,"{transcriptome}.fa"),
        profile= os.path.join(DATA_DIR,"{database}.cm.i1i")
    output:
        tblout = os.path.join(DATA_DIR,"{transcriptome}-infernal-tblout.txt"),
    log:
        "logs/{transcriptome}_cmscan.log"
    params:
        extra=config["hmmsearch"]["params"].get("extra", ""),
    threads: 4
    conda: f"file://{__path__}/wrappers/infernal/environment.yml"
    wrapper: f"file://{__path__}/wrappers/infernal/cmscan.wrapper.py"

rule run_busco:
    input: os.path.join(DATA_DIR, "{transcriptome}.fa") 
    output:
        "txome_busco/full_table_txome_busco.tsv",
    log:
        "logs/transcriptome_busco.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage_path=os.path.join(DATA_DIR, "example-busco-lineage"),
        # optional parameters
        extra=config["busco"]["params"].get("extra", ""),
    conda: f"file://{__path__}/wrappers/infernal/environment.yml"
    wrapper: f"file://{__path__}/wrappers/infernal/busco.wrapper.py"

#def parse_busco_summary:
