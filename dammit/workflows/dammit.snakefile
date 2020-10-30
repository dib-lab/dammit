# main dammit annotate snakefile
import os, sys
import numpy as np
import pandas as pd
import glob

from dammit.meta import __path__

'''
Set up paths shared by annotate and databases commands.

For databases, the results dir will be the temp folder
generated for this run; it will contain the run config.yml
and the log files.

For annotate, the results dir will be the output folder
for the annotation routine: either the default folder
based on the input transcriptome name, or the one
provided by --output-dir.
'''
results_dir = config["output_dir"]

logs = config.get("logs_dir", "logs")
logs_dir = os.path.join(results_dir, logs)

# Benchmarks aren't run for now, but we could easily turn them on
benchmarks = config.get("benchmark_dir", "benchmarks")
benchmarks_dir = os.path.join(results_dir, benchmarks)

database_dir = config['database_dir']

# we have to constrain the transcriptome name wildcard
# to keep it from slurping up intermediate files
wildcard_constraints:
    transcriptome = config['transcriptome_name'],
    database = "(?!x\.).+"


onstart:
    print("------------------------------")
    print("Just annotate it, dammit!") 
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")


onerror:
    print("Nope.\n")


'''
Somewhat nonstandard: we don't include the databases
rules for the annotation pipeline, even though the databases
outputs are inputs to the annotation rules. We want annotation
to fail if the databases haven't been installed in the provided
directory; this helps prevent users new to the software, or to
CLI bioinformatics in general, from installing multiple copies of
the databases.
'''
if config['command'] == 'databases':
    include: "databases/databases.snakefile"
else:
    include: "annotate/annotate.snakefile"
    include: "file_conversions.snakefile"
