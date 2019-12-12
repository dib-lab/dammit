# main dammit annotate snakefile
import os, sys
import numpy as np
import pandas as pd
import glob


# read in configfiles
databases = ["Pfam-A", "Rfam", "uniprot_sprot", "uniref90", "nr", "OrthoDB",  "orthodb_genes", "BUSCO"]

if quick_annotate:
    annotate_programs = ["busco", "hmmer", "infernal", "last", "transdecoder"]
    databases = ["Pfam-A", "Rfam", "uniprot_sprot", "OrthoDB"]
elif full_annotate:
    annotate_programs = ["busco", "hmmer", "infernal", "last", "transdecoder"]
    databases += ["uniref90"]
elif nr_annotate:
    annotate_programs = ["busco", "hmmer", "infernal", "last", "transdecoder"]
    databases += ["nr"]
else:
    annotate_programs = ["busco", "hmmer", "infernal", "last", "transdecoder"]
    databases = []


# generate annotation targets
output_suffixes = [config[prog]["output_suffix"] for prog in annotate_programs]
annotate_targs = [os.path.join(OUT_DIR, config[transcriptome] + suffix) for suffix in output_suffixes]

# generate database targets

# if we add the databases config into main config, can do this:
database_targs = [config["Databases"][db].get("filename","") for db in databases]

if BUSCO in databases:
    busco_dbinfo = config["Databases"]["BUSCO"] #get busco database info
    database_targs += [busco_dbinfos[db].get("filename", "") for db in busco_databases]


onstart:
    print("------------------------------")
    print("Just annotate it, dammit!") 
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")


onerror:
    print("Nope.\n")


rule all:
    input: annotation_targs

subworkflow databases:
    workdir:
        "."
    snakefile:
        "workflows/databases/databases.snakefile"
    configfile:
        "databases.json"

subworkflow annotate:
    workdir:
        "."
    snakefile:
        "workflows/databases/annotation.snakefile"
    configfile:
        "config.yml"

rule annotate:
    input:
        databases("test.txt")

#rule databases:
#    input: database_targs
