"""Snakemake wrapper for BUSCO assessment"""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2020, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
from os import path
from configparser import ConfigParser

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
config = snakemake.params.get("config", None)
if config:
    config_cmd = f" --config {config} "
else:
    config_cmd = ""
mode = snakemake.params.get("mode")
assert mode is not None, "please input a run mode: genome, transcriptome or proteins"
lineage = snakemake.params.get("lineage")
auto_lineage = snakemake.params.get("auto_lineage") # prok, euk, all
default_config = snakemake.params.get("default_configfile")
out_path= snakemake.params.get("out_path")
database_directory= snakemake.params.get("database_directory")

# modify the configfile
configur = ConfigParser()
if default_config:
    config = iniFile(default_config)
    # set path for database downloads
    if database_directory:
        configur.set("busco","download_path", os.path.abspath(out_path))
    # set path for output files
    if out_path:
        configur.set("busco_run","out_path", os.path.abspath(out_path))

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

# busco does not allow you to direct output location: handle this by moving output
outdir = str(snakemake.output[0])
if "/" in outdir:
    out_name = path.basename(outdir)
else:
    out_name = outdir

# note: --force allows snakemake to handle rewriting files as necessary
# without needing to specify *all* busco outputs as snakemake outputs
shell(
    "busco --in {snakemake.input} --out {out_name} --force "
    " --cpu {snakemake.threads} --mode {mode} {lineage_cmd} "
    " {config_cmd} {extra} {log}"
)
