import os
from dammit.meta import __path__

rule dammit_rename_transcriptome:
    input:
        config["input_transcriptome"],
    output:
        fasta=os.path.join(results_dir, "{transcriptome}.fasta"),
        names=os.path.join(results_dir, "{transcriptome}.namemap.csv")
    log:
        os.path.join(logs_dir, "{transcriptome}.rename.log")
    params:
        basename = config.get("basename", "Txome"),
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'

rule dammit_cmscan_to_gff:
    message: "Given raw input from Infernal's cmscan, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir,"{transcriptome}.{database}.cmscan-tblout.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.cmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.cmscan-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir,"{w.database}.cm.i1i")
    shell:
        """
        dammit cmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_hmmsearch_to_gff:
    message: "Given raw input from hmmer's hmmsearch, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.hmmsearch-domtbl.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.hmmsearch.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.hmmsearch-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir, '{database}.hmm.h3f') 
    shell:
        """
        dammit hmmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_hmmscan_to_gff:
    message: "Given raw input from hmmer's hmmscan, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.hmmscan-domtbl.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.hmmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.hmmscan-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir, '{database}.hmm.h3f') 
    shell:
        """
        dammit hmmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_maf_to_gff:
    message: 
        """
        Given either a raw MAF file or a CSV file the proper MAF colums,
        convert it to GFF3 and save the results.
        """
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.maf")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.lastal.maf-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir, "{database}.fasta.prj")
    shell:
        """
        dammit maf-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_shmlast_to_gff:
    message: "Given the CSV output from shmlast, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast_crbl.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast_crbl.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.shmlast-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir, "{w.database}.") #NOTE: add shmlast db extension
    shell:
        """
        dammit shmlast-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_merge_gff:
    #message: "Merge GFF files from the individual annotation programs"
    input: config["gff_files"],
    output:
        os.path.join(results_dir, "{transcriptome}.dammit.gff3"),
    log:
        os.path.join(logs_dir, "{transcriptome}.merge_gffs.log")
    shell:
        """
        dammit merge-gff3 {input} {output} 2> {log}
        """

rule dammit_annotate_fasta:
    input:
        fasta=rules.dammit_rename_transcriptome.output.fasta,
        gff=rules.dammit_merge_gff.output
    output:
        os.path.join(results_dir, "{transcriptome}.dammit.fasta")
    log:
        os.path.join(logs_dir, "{transcriptome}.annotate_fasta.log")
    shell:
        """
        dammit annotate-fasta {input.fasta} {input.gff} {output} 2> {log}
        """
