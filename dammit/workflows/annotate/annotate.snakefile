import os
from dammit.meta import __path__, __wrappers__

GLOBAL_EVALUE = config['global_evalue']
THREADS_PER_TASK = config['max_threads_per_task']

rule transdecoder_longorfs:
    message: "Run TransDecoder.LongOrfs, which fings the longest likely open reading frames."
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta')
    output:
        os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep')
    log:
        os.path.join(logs_dir, '{transcriptome}.transdecoder-longorfs.log')
    params:
        extra = config['transdecoder_longorfs']['params'].get('extra', '-m 80 ')
    threads: 1
    wrapper:
        f'file://{__wrappers__}/transdecoder/transdecoder-longorfs.wrapper.py'


# probably want to switch to hmmsearch instead of hmmscan
rule hmmscan:
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        profile = os.path.join(db_dir, 'Pfam-A.hmm.h3f'),
    output:
        # only one of these is required
        domtblout  = os.path.join(results_dir,'{transcriptome}.x.Pfam-A.hmmscan-domtbl.txt') # save parseable table of per-domain hits to file <f>
    log:
        os.path.join(logs_dir, '{transcriptome}.x.Pfam-A.hmmscan.log')
    params:
        evalue_threshold = GLOBAL_EVALUE if GLOBAL_EVALUE is not None else config['hmmscan']['params'].get("evalue", 0.00001),
        extra = config['hmmscan']['params'].get('extra', ''),
    threads: THREADS_PER_TASK
    wrapper:
        f'file://{__wrappers__}/hmmer/hmmscan.wrapper.py'


rule transdecoder_predict:
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta'),
        longorfs = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        pfam_hits = os.path.join(results_dir, '{transcriptome}.x.Pfam-A.hmmscan-domtbl.txt')
    output:
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.bed'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.cds'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.pep'),
        os.path.join(results_dir, '{transcriptome}.fasta.transdecoder.gff3')
    log:
        os.path.join(logs_dir, '{transcriptome}.transdecoder-predict.log')
    params:
        extra= config['transdecoder_predict']['params'].get('extra', '')
    threads: 1
    wrapper:
        f'file://{__wrappers__}/transdecoder/transdecoder-predict.wrapper.py'


rule hmmer_remap:
    input:
        tbl = os.path.join(results_dir,'{transcriptome}.x.Pfam-A.hmmscan-domtbl.txt'),
        pep = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep')
    output:
        os.path.join(results_dir, '{transcriptome}.x.{database}.remapped.csv')
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.hmmer_remap.log')
    threads: 1
    shell:
        """
        dammit remap-hmmer-coords {input.tbl} {input.pep} {output} 2> {log}
        """


rule lastal:
    input:
        data = os.path.join(results_dir, '{transcriptome}.fasta'),
        lastdb = os.path.join(db_dir, '{database}.fasta.prj')
    output:
        maf = os.path.join(results_dir, '{transcriptome}.x.{database}.lastal.maf')
    params:
        evalue_threshold = GLOBAL_EVALUE if GLOBAL_EVALUE is not None else config['lastal']['params'].get('evalue', 0.00001),
        frameshift_cost = config['lastal']['params'].get('frameshift_cost', 15),
        extra           = config['lastal']['params'].get('extra', ''),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.lastal.log')
    threads: THREADS_PER_TASK
    wrapper: 
        f'file://{__wrappers__}/last/lastal.wrapper.py'


rule shmlast_crbl:
    input:
        query = os.path.join(results_dir, '{transcriptome}.fasta'),
        database = lambda w: config["user_dbs"][w.database] # get full path from dictionary in configfile 
    output:
        os.path.join(results_dir, '{transcriptome}.x.{database}.shmlast_crbl.csv')
    params:
        search_type="crbl",
        evalue = GLOBAL_EVALUE if GLOBAL_EVALUE is not None else config['shmlast']['params'].get('evalue', 0.00001),
        extra = config['shmlast']['params'].get('extra', ''),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.shmlast.log')
    threads: THREADS_PER_TASK
    wrapper: f'file://{__wrappers__}/shmlast/shmlast.wrapper.py'


rule cmscan:
    input:
        fasta   = os.path.join(results_dir,'{transcriptome}.fasta'),
        profile = os.path.join(db_dir,'{database}.cm.i1i')
    output:
        tblout = os.path.join(results_dir,'{transcriptome}.x.{database}.cmscan-tblout.txt'),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.cmscan.log')
    params:
        evalue_threshold = GLOBAL_EVALUE if GLOBAL_EVALUE is not None else config['cmscan']['params'].get('evalue', 0.00001),
        extra = config['hmmsearch']['params'].get('extra', ''),
    threads: THREADS_PER_TASK
    wrapper: f'file://{__wrappers__}/infernal/cmscan.wrapper.py'


rule busco_transcripts:
    input:
        fasta = os.path.join(results_dir,'{transcriptome}.fasta')
    output:
        os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', 'run_{busco_db}', 'short_summary.txt')
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.log")
    benchmark:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.benchmark")
    threads: THREADS_PER_TASK
    params:
        out_name = '{busco_db}_outputs',
        out_path = lambda w: os.path.join(results_dir, f'{w.transcriptome}.busco'),
        config = config['busco']['configfile'],
        mode = "transcriptome",
        lineage = lambda w: w.busco_db,
        database_directory = db_dir,
        #auto_lineage='euk', # enabled in wrapper, but not using this bc it changes output dir structure
        extra = config['busco']['params'].get('extra', ''),
    wrapper: f'file://{__wrappers__}/busco/busco.wrapper.py'


localrules: plot_busco_summaries


def expand_busco_summaries(w):
    return expand(os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', "run_{busco_db}", "short_summary.txt"), 
                         busco_db = config["busco_groups"],
                         transcriptome = w.transcriptome)


rule plot_busco_summaries:
    input: expand_busco_summaries
    output: os.path.join(results_dir, '{transcriptome}.busco', "summary_figure.png")
    log: os.path.join(logs_dir, "{transcriptome}.busco", "summary_figure.log")
    conda:
        f'file://{__path__}/wrappers/busco/environment.yaml'
    params:
        summary_dir = lambda w: os.path.join(results_dir, f'{w.transcriptome}.busco', 'summary_data')
    threads: 1
    shell:
        """
        mkdir -p {params.summary_dir}
        cp {wildcards.transcriptome}.busco/*_outputs/short_summary.*.txt {params.summary_dir}
        generate_plot.py --working_directory {params.summary_dir} 2> {log}
        cp {params.summary_dir}/busco_figure.png {output}
        """
