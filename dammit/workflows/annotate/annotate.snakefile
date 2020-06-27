import os
from dammit.meta import __path__, __wrappers__

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
    threads: 4
    wrapper:
        f'file://{__wrappers__}/transdecoder/transdecoder-longorfs.wrapper.py'


'''
rule hmmsearch:
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        profile = os.path.join(db_dir, 'Pfam-A.hmm.h3f')
    output:
        # only one of these is required
        tblout = os.path.join(results_dir, '{transcriptome}.x.Pfam-A.hmmsearch-domtbl.txt'), # save parseable table of per-domain hits to file <f>
        #tblout=os.path.join(results_dir, '{transcriptome}_{database}.hmmsearch-tbl.txt'), # save parseable table of per-sequence hits to file <f>
        #alignment_hits=os.path.join(results_dir,'{transcriptome}_{database}.hmmsearch-alignment-hits.txt'), # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
        #outfile=os.path.join(results_dir,'{transcriptome}_{database}.hmmsearch-out.txt'), # Direct the main human-readable output to a file <f> instead of the default stdout. 
    log:
        os.path.join(logs_dir, '{transcriptome}.x.Pfam-A.hmmsearch.log')
    params:
        evalue_threshold = config['hmmscan']['params'].get("evalue", 0.00001),
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra = config['hmmsearch']['params'].get('extra', ''),
    threads: 4
    wrapper:
        f'file://{__wrappers__}/hmmer/hmmsearch.wrapper.py'
'''


# probably want to switch to hmmsearch instead of hmmscan
rule hmmscan:
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        profile = os.path.join(db_dir, 'Pfam-A.hmm.h3f'),
    output:
        # only one of these is required
        domtblout  = os.path.join(results_dir,'{transcriptome}.x.Pfam-A.hmmscan-domtbl.txt'), # save parseable table of per-domain hits to file <f>
        outfile    = os.path.join(results_dir,'{transcriptome}.x.Pfam-A.hmmscan-out.txt'), # Direct the main human-readable output to a file <f> instead of the default stdout.
    log:
        os.path.join(logs_dir, '{transcriptome}.x.Pfam-A.hmmscan.log')
    params:
        evalue_threshold = config['hmmscan']['params'].get("evalue", 0.00001),
        # if bitscore threshold provided, hmmscan will use that instead
        #score_threshold=50,
        extra = config['hmmscan']['params'].get('extra', ''),
    threads: 4
    wrapper:
        f'file://{__wrappers__}/hmmer/hmmscan.wrapper.py'


rule transdecoder_predict:
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta'),
        longorfs = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        pfam = os.path.join(results_dir, '{transcriptome}.x.Pfam-A.hmmscan-domtbl.txt')
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
    wrapper:
        f'file://{__wrappers__}/transdecoder/transdecoder-predict.wrapper.py'


rule lastal:
    input:
        data = os.path.join(results_dir, '{transcriptome}.fasta'),
        lastdb = os.path.join(db_dir, '{database}.fasta.prj')
    output:
        maf = os.path.join(results_dir, '{transcriptome}.x.{database}.lastal.maf')
    params:
        frameshift_cost = config['lastal']['params'].get('frameshift_cost', 15),
        extra           = config['lastal']['params'].get('extra', ''),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.lastal.log')
    threads: 8
    wrapper: f'file://{__wrappers__}/last/lastal.wrapper.py'

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
        extra = config['hmmsearch']['params'].get('extra', ''),
    threads: 4
    wrapper: f'file://{__wrappers__}/infernal/cmscan.wrapper.py'

rule busco_transcripts:
    input:
        fasta=os.path.join(results_dir,'{transcriptome}.fasta'),
        config=config["busco_config_file"]
    output:
        os.path.join(results_dir, '{transcriptome}.busco.{busco_db}', "run_{busco_db}", "short_summary.txt")
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.log")
    benchmark:
        os.path.join(logs_dir, "{transcriptome}.x.{busco_db}.benchmark")
    threads: 8
    params:
        mode = "transcriptome",
        lineage=lambda w: w.busco_db,
        database_directory= db_dir,
        #auto_lineage='euk', # enabled in wrapper, but not using this bc it changes output dir structure
        extra = config['busco']['params'].get('extra', ''),
    wrapper: f'file://{__wrappers__}/busco/busco.wrapper.py'

def aggregate_busco_summaries(w):
    busco_files = expand(os.path.join(results_dir, '{transcriptome}.busco.{busco_db}', "run_{busco_db}", "short_summary.txt"), busco_db = config["busco_groups"], transcriptome=w.transcriptome)
    summary_files=[]
    for s in busco_files:
        buscoD = os.path.dirname(os.path.dirname(s))
        summary_files+= glob.glob(os.path.join(buscoD, "short_summary*.txt"))
    #  (short_summary.[generic|specific].dataset.label.txt)
    return summary_files

localrules: plot_busco_summaries

rule plot_busco_summaries:
    input: aggregate_busco_summaries
    output: os.path.join(results_dir, "{transcriptome}.busco_summary_plot.png")
    log: os.path.join(logs_dir, "busco", "{transcriptome}.x.plot_summaries.log")
    params:
        outdir= os.path.join(results_dir, "busco_summaries")
    conda:
        f'file://{__path__}/wrappers/busco/environment.yaml'
    shell:
        """
        mkdir -p {params.outdir}
        cp {input} {params.outdir}
        generate_plot.py --working_directory {params.outdir} 2> {log}
        cp {params.outdir}/busco_figure.png {output}
        """
