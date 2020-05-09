import os
from dammit.meta import __path__

rule transdecoder_longorfs:
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta')
    output:
        os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep')
    log:
        os.path.join(logs_dir, '{transcriptome}.transdecoder-longorfs.log')
    params:
        extra = config['transdecoder_longorfs']['params'].get('extra', '-m 80 ')
    threads: 4
    conda: 
        f'file://{__path__}/wrappers/transdecoder/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/transdecoder/transdecoder-longorfs.wrapper.py'

rule transdecoder_predict:
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta'),
        longorfs = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep')
    output:
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.bed'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.cds'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.pep'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.gff3')
    log:
        os.path.join(logs_dir, '{transcriptome}.transdecoder-predict.log')
    params:
        extra= config['transdecoder_predict']['params'].get('extra', '')
    threads: 4
    conda: 
        f'file://{__path__}/wrappers/transdecoder/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/transdecoder/transdecoder-predict.wrapper.py'

rule lastal:
    input:
        data = os.path.join(results_dir, '{transcriptome}.fasta'),
        lastdb = os.path.join(db_dir, '{database}.fasta.prj')
    output:
        maf = os.path.join(results_dir, '{transcriptome}.{database}.lastal.maf')
    params:
        frameshift_cost = config['lastal']['params'].get('frameshift_cost', 15),
        extra           = config['lastal']['params'].get('extra', ''),
    log:
        os.path.join(logs_dir, '{transcriptome}.{database}.lastal.log')
    threads: 8
    conda: f'file://{__path__}/wrappers/last/environment.yml'
    script: f'file://{__path__}/wrappers/last/lastal.wrapper.py'

rule shmlast_crbl:
    input:
        query = os.path.join(results_dir, '{transcriptome}.fasta'),
        database = lambda w: config["user_dbs"][w.database] # get full path from dictionary in configfile 
    output:
        os.path.join(results_dir, '{transcriptome}.x.{database}.shmlast_crbl.csv')
    params:
        search_type="crbl",
        evalue = config['shmlast']['params'].get('evalue', ""),
        extra = config['shmlast']['params'].get('extra', ''),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.shmlast.log')
    threads: 8
    conda: f'file://{__path__}/wrappers/shmlast/environment.yml'
    script: f'file://{__path__}/wrappers/shmlast/shmlast.wrapper.py'

# probably want to switch to hmmsearch instead of hmmscan
rule hmmscan:
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.fasta'),
        profile = os.path.join(db_dir, '{database}.hmm.h3f'),
    output:
        # only one of these is required
        tblout     = os.path.join(results_dir, '{transcriptome}.hmmscan-tbl.txt'), # save parseable table of per-sequence hits to file <f>
        domtblout  = os.path.join(results_dir,'{transcriptome}.hmmscan-domtbl.txt'), # save parseable table of per-domain hits to file <f>
        pfamtblout = os.path.join(results_dir,'{transcriptome}.hmmscan-pfamtbl.txt'), # save table of hits and domains to file, in Pfam format <f>
        outfile    = os.path.join(results_dir,'{transcriptome}.hmmscan-out.txt'), # Direct the main human-readable output to a file <f> instead of the default stdout.
    log:
        os.path.join(logs_dir, '{transcriptome}_hmmscan.log')
    params:
        evalue_threshold = config['hmmscan']['params'].get("evalue", 0.00001),
        # if bitscore threshold provided, hmmscan will use that instead
        #score_threshold=50,
        extra = config['hmmscan']['params'].get('extra', ''),
    threads: 4
    conda:
        f'file://{__path__}/wrappers/hmmer/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/hmmer/hmmscan.wrapper.py'

rule hmmsearch:
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.fasta'),
        profile = os.path.join(db_dir, '{database}.hmm.h3f')
    output:
        # only one of these is required
        domtblout = os.path.join(results_dir, '{transcriptome}_{database}.hmmsearch-domtbl.txt'), # save parseable table of per-domain hits to file <f>
        #tblout=os.path.join(results_dir, '{transcriptome}_{database}.hmmsearch-tbl.txt'), # save parseable table of per-sequence hits to file <f>
        #alignment_hits=os.path.join(results_dir,'{transcriptome}_{database}.hmmsearch-alignment-hits.txt'), # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
        #outfile=os.path.join(results_dir,'{transcriptome}_{database}.hmmsearch-out.txt'), # Direct the main human-readable output to a file <f> instead of the default stdout. 
    log:
        os.path.join(logs_dir, '{transcriptome}_{database}_hmmsearch.log')
    params:
        evalue_threshold = config['hmmscan']['params'].get("evalue", 0.00001),
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra = config['hmmsearch']['params'].get('extra', ''),
    threads: 4
    conda:
        f'file://{__path__}/wrappers/hmmer/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/hmmer/hmmsearch.wrapper.py'

rule cmscan:
    input:
        fasta   = os.path.join(results_dir,'{transcriptome}.fasta'),
        profile = os.path.join(db_dir,'{database}.cm.i1i')
    output:
        tblout = os.path.join(results_dir,'{transcriptome}.{database}.cmscan-tblout.txt'),
    log:
        os.path.join(logs_dir, '{transcriptome}.{database}.cmscan.log')
    params:
        extra = config['hmmsearch']['params'].get('extra', ''),
    threads: 4
    conda:
        f'file://{__path__}/wrappers/infernal/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/infernal/cmscan.wrapper.py'

rule busco_transcripts:
    input: os.path.join(results_dir,'{transcriptome}.fasta'),
    output:
        directory(os.path.join(results_dir, '{transcriptome}.busco.{busco_db}')),
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.log")
    benchmark:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.benchmark")
    threads: 8
    params:
        mode = "transcriptome",
        lineage=lambda w: w.busco_db,
        config=config.get("busco_config_file", None), # enable user to input custom config
        default_config = os.path.join(db_dir, '{database}.ini'), # or hardcode busco_config.ini?
        database_directory= db_dir,
        #auto_lineage='euk', # enabled in wrapper, but not using this bc it changes output dir structure
        extra = config['busco']['params'].get('extra', ''),
    conda:
        f'file://{__path__}/wrappers/busco/environment.yml'
    wrapper:
        f'file://{__path__}/wrappers/busco/busco.wrapper.py'
