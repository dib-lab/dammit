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
    params:
        database=lambda w: os.path.join(db_dir,"{w.database}.cm.i1i")
    shell:
        """
        dammit convert cmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
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
        dammit convert hmmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
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
        dammit convert hmmscan-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
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
        dammit convert maf-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """

rule dammit_shmlast_to_gff:
    message: "Given the CSV output from shmlast, convert it to GFF3 and save the results."
    input: 
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast.csv")
    output:
        os.path.join(results_dir, "{transcriptome}.{database}.shmlast.gff3")
    log:
        os.path.join(logs_dir, "{transcriptome}.{database}.shmlast-to-gff3.log")
    params:
        database=lambda w: os.path.join(db_dir, "{w.database}.") #NOTE: add shmlast db extension
    shell:
        """
        dammit convert shmlast-to-gff3 --dbxref {params.database} {input} {output} 2> {log}
        """
