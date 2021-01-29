import os
from dammit.meta import __path__, __wrappers__

GLOBAL_EVALUE = config['global_evalue']
THREADS_PER_TASK = config['max_threads_per_task']


rule dammit_rename_transcriptome:
    message:
        """
        Reformat and rename FASTA headers for compatibility. 
        """
    input:
        config["input_transcriptome"],
    output:
        fasta=os.path.join(results_dir, "{transcriptome}.fasta"),
        names=os.path.join(results_dir, "{transcriptome}.namemap.csv")
    log:
        os.path.join(logs_dir, "{transcriptome}.rename.log")
    threads: 1
    params:
        basename = config.get("basename", "Txome"),
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'


rule dammit_transcriptome_stats:
    message:
        """
        Calculate some basic stats and metrics on the input transcriptome.
        """
    input:
        fasta = os.path.join(results_dir, '{transcriptome}.fasta')
    output:
        stats_fn = os.path.join(results_dir, '{transcriptome}.stats.json'),
        per_transcript_fn = os.path.join(results_dir, '{transcriptome}.per-transcript-stats.csv')
    log:
        os.path.join(logs_dir, '{transcriptome}.stats.log')
    threads: 1
    shell:
        """
        dammit transcriptome-stats {input.fasta} {output.stats_fn} {output.per_transcript_fn}
        """


rule transdecoder_longorfs:
    message: 
        """
        Run TransDecoder.LongOrfs, which finds the longest likely open reading frames.
        """
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
    message:
        """
        Run hmmscan against Pfam-A on the ORFs produced by TransDecoder.LongOrfs.
        """
    input:
        fasta   = os.path.join(results_dir, '{transcriptome}.transdecoder_dir/longest_orfs.pep'),
        profile = os.path.join(database_dir, 'Pfam-A.hmm.h3f'),
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
    message:
        """
        Run TransDecoder.Predict, using the ORFs and the Pfam-A domains, to predict transcript features.
        """
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
    message:
        """
        Remap the coordinates from the ORF-to-Pfam-A domain predictions back to their
        locations on the input transcripts.
        """
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
    message:
        """
        Find best-hits between the transcripts and the given protein database.
        """
    input:
        data = os.path.join(results_dir, '{transcriptome}.fasta'),
        lastdb = os.path.join(database_dir, '{database}.fasta.prj')
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
    message:
        """
        Find "conditional reciprocal best hits" between the transcriptome and the
        user-provided protein databases. See the docs for details on this method.
        """
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
    message:
        """
        Use Infernal's cmscan to search for non-coding RNAs using the
        Rfam secondary-structure covariance model database.
        """
    input:
        fasta   = os.path.join(results_dir,'{transcriptome}.fasta'),
        profile = os.path.join(database_dir,'{database}.cm.i1i')
    output:
        tblout = os.path.join(results_dir,'{transcriptome}.x.{database}.cmscan-tblout.txt'),
    log:
        os.path.join(logs_dir, '{transcriptome}.x.{database}.cmscan.log')
    params:
        evalue_threshold = GLOBAL_EVALUE if GLOBAL_EVALUE is not None else config['cmscan']['params'].get('evalue', 0.00001),
        extra = config['cmscan']['params'].get('extra', ''),
    threads: THREADS_PER_TASK
    wrapper: f'file://{__wrappers__}/infernal/cmscan.wrapper.py'


rule busco_transcripts:
    message:
        """
        Run BUSCO to assess the completeness of the transcriptome assembly
        using a set of benchmarking single-copy orthologs.
        """
    input:
        fasta = os.path.join(results_dir,'{transcriptome}.fasta')
    output:
        os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', 'run_{busco_db}', 'short_summary.txt'),
        os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', 'run_{busco_db}', 'full_table.tsv')
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
        database_directory = database_dir,
        #auto_lineage='euk', # enabled in wrapper, but not using this bc it changes output dir structure
        extra = config['busco']['params'].get('extra', ''),
    wrapper: f'file://{__wrappers__}/busco/busco.wrapper.py'


localrules: plot_busco_summaries


def expand_busco_files(w, basename='short_summary.txt'):
    return expand(os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', "run_{busco_db}", basename), 
                         busco_db = config["busco_groups"],
                         transcriptome = w.transcriptome)


rule plot_busco_summaries:
    message:
        """
        Plot the BUSCO results for the user-provided lineage databases.
        """
    input: expand_busco_files
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


rule dammit_busco_to_gff:
    message:
        """
        Convert BUSCO results to a GFF3 representation.
        """
    input: 
        busco_fn = os.path.join(results_dir, '{transcriptome}.busco', '{busco_db}_outputs', 'run_{busco_db}', 'full_table.tsv'),
        lens_fn = os.path.join(results_dir, '{transcriptome}.per-transcript-stats.csv')
    output:
        os.path.join(results_dir, '{transcriptome}.x.busco.{busco_db}.gff3')
    log:
        os.path.join(logs_dir, '{transcriptome}.{busco_db}.busco-to-gff3.log')
    threads: 1
    shell:
        """
        dammit busco-to-gff3 {input.busco_fn} {input.lens_fn} {output}
        """


rule dammit_cmscan_to_gff:
    message: 
        """
        Given raw input from Infernal's cmscan, convert it to GFF3 and save the results.
        """
    input: 
        os.path.join(results_dir,"{transcriptome}.x.{database}.cmscan-tblout.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.x.{database}.cmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{database}.cmscan-to-gff3.log")
    threads: 1
    shell:
        """
        dammit cmscan-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """


rule dammit_hmmscan_to_gff:
    message:
        """
        Convert hmmscan's domain table output to GFF3.
        """
    input: 
        os.path.join(results_dir, "{transcriptome}.x.{database}.remapped.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.x.{database}.hmmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{database}.hmmscan-to-gff3.log")
    threads: 1
    shell:
        """
        dammit hmmscan-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """


rule dammit_maf_best_hits:
    message:
        """
        Filter out only the best (top-scoring) hits for each transcript from the 
        MAF alignments.
        """
    input:
        os.path.join(results_dir, "{transcriptome}.x.{database}.lastal.maf")
    output:
        os.path.join(results_dir, "{transcriptome}.x.{database}.lastal.best.csv")
    threads: 1
    shell:
        """
        dammit best-hits {input} {output}
        """


rule dammit_maf_to_gff:
    message: 
        """
        Given either a raw MAF file or a CSV file the proper MAF colums,
        convert it to GFF3 and save the results.
        """
    input: 
        os.path.join(results_dir, "{transcriptome}.x.{database}.lastal.best.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.x.{database}.lastal.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{database}.lastal.maf-to-gff3.log")
    threads: 1
    shell:
        """
        dammit maf-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """

rule dammit_shmlast_to_gff:
    message: 
        """
        Given the CSV output from shmlast, convert it to GFF3 and save the results.
        """
    input: 
        os.path.join(results_dir, "{transcriptome}.x.{database}.shmlast_crbl.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.x.{database}.shmlast_crbl.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.x.{database}.shmlast-to-gff3.log")
    threads: 1
    shell:
        """
        dammit shmlast-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """

rule dammit_merge_gff:
    message: 
        """
        Merge GFF files from the individual annotation programs.
        """
    input: config["gff_files"],
    output:
        os.path.join(results_dir, "{transcriptome}.dammit.gff3"),
    log:
        os.path.join(logs_dir, "{transcriptome}.merge_gffs.log")
    threads: 1
    shell:
        """
        dammit merge-gff3 {input} {output} 2> {log}
        """

rule dammit_annotate_fasta:
    message:
        """
        Annotate the headers of a FASTA file with a summary of each sequence.
        """
    input:
        fasta=rules.dammit_rename_transcriptome.output.fasta,
        gff=rules.dammit_merge_gff.output
    output:
        os.path.join(results_dir, "{transcriptome}.dammit.fasta")
    log:
        os.path.join(logs_dir, "{transcriptome}.annotate_fasta.log")
    threads: 1
    shell:
        """
        dammit annotate-fasta {input.fasta} {input.gff} {output} 2> {log}
        """
