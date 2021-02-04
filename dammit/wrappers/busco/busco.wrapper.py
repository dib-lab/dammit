"""Snakemake wrapper for BUSCO assessment"""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2020, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell
from configparser import ConfigParser

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
mode = snakemake.params.get("mode")
assert mode is not None, "please input a run mode: genome, transcriptome or proteins"

evalue = snakemake.params.get('evalue', None)
lineage = snakemake.params.get("lineage")
auto_lineage = snakemake.params.get("auto_lineage") # prok, euk, all
database_directory = snakemake.params.get("database_directory")
config = snakemake.params.get("config", None)


out_path = snakemake.params.get("out_path", None)
assert out_path is not None, "must specific out_path param for busco"
out_name = snakemake.params.get("out_name", None)
assert out_name is not None, "must specify an out_name param for busco"
assert snakemake.output[0].startswith(out_path)

# handle config file
config_cmd = ""
if config and database_directory:
    configur = ConfigParser()
    config = configur.read(config)
    # set path for database downloads
    if database_directory:
        configur.set("busco_run","download_path", database_directory)

    # print configfile to output directory
    configfile = os.path.join(out_path, "run.busco_config.ini")
    with open(configfile, "w") as outF:
        configur.write(outF)
    # cmd to point busco to this new configfile
    config_cmd = f" --config {configfile} "

#assert lineage is not None, "please input the path to a lineage for busco assessment"
if lineage is not None:
    lineage_cmd = f" -l {lineage} "
elif auto_lineage is not None:
    if auto_lineage == "prok":
        lineage_cmd = " --auto-lineage-prok "
    elif auto_lineage == "euk":
        lineage_cmd = " --auto-lineage-euk "
else:   # doesn't matter if auto-lineage is all or left blank. default to auto if nothing else is provided
    lineage_cmd = " --auto-lineage "

evalue_cmd = '--evalue ' + str(evalue) if evalue is not None else ''

# note: --force allows snakemake to handle rewriting files as necessary
# without needing to specify *all* busco outputs as snakemake outputs
shell(
    "busco --in {snakemake.input.fasta} --out_path {out_path} --out {out_name} --force "
    " --cpu {snakemake.threads} {evalue_cmd} --mode {mode} {lineage_cmd} "
    " {config_cmd} {extra} {log}"
)