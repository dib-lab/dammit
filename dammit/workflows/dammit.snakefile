'''
    | |                         (_) |  
  __| | __ _ _ __ ___  _ __ ___  _| |_ 
 / _` |/ _` | '_ ` _ \| '_ ` _ \| | __|
| (_| | (_| | | | | | | | | | | | | |_ 
 \__,_|\__,_|_| |_| |_|_| |_| |_|_|\__|

(c) Camille Scott 2015-2020, N. Tessa Pierce 2018-2020
MIT License 

This is the dammit Snakemake entry point. This workflow is meant to be
executed through dammit, not as a standalone Snakefile; advanced users could do
this, but they would need to provide their own config file and supply all the
metadata directories correctly.

This workflow takes care of the databases installation and the annotation
pipeline in separate Snakefiles. They are meant to be run independently.
'''

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


onsuccess:
    print("\n🟢 Workflow completed!")


onerror:
    print("\n🔴 Workflow encountered an error.")


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
    include: "databases.snakefile"
else:
    include: "annotate.snakefile"
