# main dammit annotate snakefile
import os, sys
import numpy as np
import pandas as pd
import glob

from dammit.meta import __path__

results_dir = config["dammit_dir"]
logs = config.get("logs_dir", "logs")
logs_dir = os.path.join(results_dir, logs)

benchmarks = config.get("benchmark_dir", "benchmarks")
benchmarks_dir = os.path.join(results_dir, benchmarks)
db_dir = config['db_dir']


onstart:
    print("------------------------------")
    print("Just annotate it, dammit!") 
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")


onerror:
    print("Nope.\n")


include: "setup.snakefile"
include: "file_conversions.snakefile"
include: "databases/databases.snakefile"
include: "annotate/annotate.snakefile"
#subworkflow databases:
#    workdir:
#        "."
#    snakefile:
#        "databases/databases.snakefile"
#    configfile: f"{__path__}/config.yml" 
#        "databases.yml", 
        #"config.yml"

#subworkflow annotate:
#    workdir:
#        "."
#    snakefile:
#        "annotate/annotate.snakefile"
#    configfile:
#        "config.yml"

#rule annotate:
#    input:
#        databases("test.txt")

#rule databases:
#    input: database_targs
