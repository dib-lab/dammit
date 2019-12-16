import os
from dammit.meta import __path__


# this gets replaced with dammit's rename_transcriptome
rule cp_transcriptome:
    input: config["input_transcriptome"]
    output: os.path.join(config["dammit_dir"], "{transcriptome}.fasta")
    log: os.path.join(config["dammit_dir"], 'logs/cp-{transcriptome}.log')
    shell:
        """
        cp {input} {output} 2> {log}
        """

rule transdecoder_longorfs:
    input:
        fasta = os.path.join(config["dammit_dir"], "{transcriptome}.fasta")
    output:
        os.path.join(config["dammit_dir"], "{transcriptome}.transdecoder_dir/longest_orfs.pep")
    log:
        os.path.join(config["dammit_dir"], 'logs/{transcriptome}.transdecoder-longorfs.log')
    params:
        extra=config["transdecoder_longorfs"]["params"].get("extra", "-m 80 ")
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-longorfs.wrapper.py"

rule transdecoder_predict:
    input:
        fasta = os.path.join(config["dammit_dir"], "{transcriptome}.fasta"),
        longorfs = os.path.join(config["dammit_dir"], "{transcriptome}.transdecoder_dir/longest_orfs.pep")
    output:
        os.path.join(config["dammit_dir"], "{transcriptome}.fasta.transdecoder.bed"),
        os.path.join(config["dammit_dir"], "{transcriptome}.fasta.transdecoder.cds"),
        os.path.join(config["dammit_dir"], "{transcriptome}.fasta.transdecoder.pep"),
        os.path.join(config["dammit_dir"], "{transcriptome}.fasta.transdecoder.gff3")
    log:
        os.path.join(config["dammit_dir"], 'logs/{transcriptome}.transdecoder-predict.log')
    params:
        extra= config["transdecoder_predict"]["params"].get("extra", "")
    threads: 4
    conda: 
        f"file://{__path__}/wrappers/transdecoder/environment.yml"
    wrapper:
        f"file://{__path__}/wrappers/transdecoder/transdecoder-predict.wrapper.py"

rule lastal:
    input:
        data=os.path.join(config["dammit_dir"], "{transcriptome}.fasta"),
        lastdb =os.path.join(config["db_dir"], "{database}.fa.prj")
    output:
        maf = os.path.join(config["dammit_dir"], "{transcriptome}_{database}.lastal.maf")
    params:
        frameshift_cost = config["lastal"]["params"].get("frameshift_cost", 15),
        extra=config["lastal"]["params"].get("extra", ""),
    log:
        os.path.join(config["dammit_dir"], "logs/{transcriptome}_{database}.lastal.log")
    threads: 8
    conda: f"file://{__path__}/wrappers/last/environment.yml"
    script: f"file://{__path__}/wrappers/last/lastal.wrapper.py"

# probably want to switch to hmmsearch instead of hmmscan
rule hmmscan:
    input:
        fasta=os.path.join(config["dammit_dir"], "{transcriptome}.fasta"),
        profile=os.path.join(config["db_dir"], "{database}.hmm.h3f"),
    output:
        # only one of these is required
        tblout=os.path.join(config["dammit_dir"], "{transcriptome}.hmmscan-tbl.txt"), # save parseable table of per-sequence hits to file <f>
        domtblout=os.path.join(config["dammit_dir"],"{transcriptome}.hmmscan-domtbl.txt"), # save parseable table of per-domain hits to file <f>
        pfamtblout=os.path.join(config["dammit_dir"],"{transcriptome}.hmmscan-pfamtbl.txt"), # save table of hits and domains to file, in Pfam format <f>
        outfile=os.path.join(config["dammit_dir"],"{transcriptome}.hmmscan-out.txt"), # Direct the main human-readable output to a file <f> instead of the default stdout.
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
        fasta=os.path.join(config["dammit_dir"], "{transcriptome}.fasta"),
        profile=os.path.join(config["db_dir"], "{database}.hmm.h3f")
    output:
        # only one of these is required
        domtblout=os.path.join(config["dammit_dir"],"{transcriptome}_{database}.hmmsearch-domtbl.txt"), # save parseable table of per-domain hits to file <f>
        #tblout=os.path.join(config["dammit_dir"], "{transcriptome}_{database}.hmmsearch-tbl.txt"), # save parseable table of per-sequence hits to file <f>
        #alignment_hits=os.path.join(config["dammit_dir"],"{transcriptome}_{database}.hmmsearch-alignment-hits.txt"), # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
        #outfile=os.path.join(config["dammit_dir"],"{transcriptome}_{database}.hmmsearch-out.txt"), # Direct the main human-readable output to a file <f> instead of the default stdout. 
    log:
        os.path.join(config["dammit_dir"], "logs/{transcriptome}_{database}_hmmsearch.log")
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
        fasta=os.path.join(config["dammit_dir"],"{transcriptome}.fasta"),
        profile=os.path.join(config["db_dir"],"{database}.cm.i1i")
    output:
        tblout=os.path.join(config["dammit_dir"],"{transcriptome}_{database}.cmscan-tblout.txt"),
    log:
        os.path.join(config["dammit_dir"], "logs/{transcriptome}_{database}.cmscan.log")
    params:
        extra=config["hmmsearch"]["params"].get("extra", ""),
    threads: 4
    conda: f"file://{__path__}/wrappers/infernal/environment.yml"
    wrapper: f"file://{__path__}/wrappers/infernal/cmscan.wrapper.py"

rule run_busco:
    input: 
        trans = os.path.join(config["dammit_dir"], "{transcriptome}.fa"),
        donefile = os.path.join(config["db_dir"], "{busco_db}.done")
    output:
        os.path.join(config["dammit_dir"], "txome_busco/{transcriptome}_{busco_db}_full_table.tsv"),
    log:
        "logs/{transcriptome}_{busco_db}.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage_path=os.path.join(config["db_dir"], "{busco_db}"),
        # optional parameters
        extra=config["busco"]["params"].get("extra", ""),
    conda: f"file://{__path__}/wrappers/busco/environment.yml"
    wrapper: f"file://{__path__}/wrappers/busco/busco.wrapper.py"

#def parse_busco_summary:
