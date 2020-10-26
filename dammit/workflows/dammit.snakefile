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

wildcard_constraints:
    transcriptome = config['transcriptome_name'],
    database = "(?!x\.).+"
    #database = "(?<!x\.).+"


onstart:
    print("------------------------------")
    print("Just annotate it, dammit!") 
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")


onerror:
    print("Nope.\n")


include: "databases/databases.snakefile"
include: "annotate/annotate.snakefile"
include: "file_conversions.snakefile"
