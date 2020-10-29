# main dammit annotate snakefile
import os, sys
import numpy as np
import pandas as pd
import glob

from dammit.meta import __path__


results_dir = config["output_dir"]
logs = config.get("logs_dir", "logs")
logs_dir = os.path.join(results_dir, logs)

benchmarks = config.get("benchmark_dir", "benchmarks")
benchmarks_dir = os.path.join(results_dir, benchmarks)
database_dir = config['database_dir']

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

if config['command'] == 'databases':
    include: "databases/databases.snakefile"
else:
    include: "annotate/annotate.snakefile"
    include: "file_conversions.snakefile"
