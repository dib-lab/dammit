import os
from dammit.meta import __path__

rule dammit_cmscan_to_gff:
    message: "Given raw input from Infernal's cmscan, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir,"{transcriptome}.{database}.cmscan-tblout.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.cmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.cmscan-to-gff3.log")
    shell:
        """
        dammit cmscan-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """

rule dammit_hmmsearch_to_gff:
    message: "Given raw input from hmmer's hmmsearch, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.hmmsearch-domtbl.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.hmmsearch.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.hmmsearch-to-gff3.log")
    shell:
        """
        dammit hmmscan-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """

rule dammit_hmmscan_to_gff:
    message: "Given raw input from hmmer's hmmscan, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.hmmscan-domtbl.txt")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.hmmscan.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.hmmscan-to-gff3.log")
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
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.maf")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.best.csv")
    shell:
        "dammit best-hits {input} {output}"


rule dammit_maf_to_gff:
    message: 
        """
        Given either a raw MAF file or a CSV file the proper MAF colums,
        convert it to GFF3 and save the results.
        """
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.best.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.lastal.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.lastal.maf-to-gff3.log")
    shell:
        """
        dammit maf-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """

rule dammit_shmlast_to_gff:
    message: "Given the CSV output from shmlast, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast_crbl.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast_crbl.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.shmlast-to-gff3.log")
    shell:
        """
        dammit shmlast-to-gff3 --dbxref {wildcards.database} {input} {output} 2> {log}
        """
