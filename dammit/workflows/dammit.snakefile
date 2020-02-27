# main dammit annotate snakefile
import os, sys
import numpy as np
import pandas as pd
import glob

from dammit.meta import __path__


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
