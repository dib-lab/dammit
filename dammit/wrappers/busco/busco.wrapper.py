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
mode = snakemake.params.get("mode")
assert mode is not None, "please input a run mode: genome, transcriptome or proteins"
lineage = snakemake.params.get("lineage")
auto_lineage = snakemake.params.get("auto_lineage") # prok, euk, all
database_directory= snakemake.params.get("database_directory")

config = snakemake.params.get("config", None) # user's custom config (not modified)
default_config = snakemake.params.get("default_configfile", None) # default busco config that we will modify below
#maybe add an assertion to check they're not inputting both?

# separate output directory and output name from snakemake.output
outdir = str(snakemake.output[0])
out_path = ""
if "/" in outdir:
    out_path = path.dirname(outdir)
    out_name = path.basename(outdir)
else:
    out_name = outdir

# handle config files
config_cmd = ""
if config: # enable user to put in custom config
    config_cmd = f" --config {config} "
elif default_config: # modify with output path and database location
    configur = ConfigParser()
    config = iniFile(default_config)
    # set path for database downloads
    if database_directory:
        configur.set("busco","download_path", os.path.abspath(out_path))
    # set path for output files
    if out_path:
        configur.set("busco_run","out_path", os.path.abspath(out_path))

    #print configfile to output directory
    configfile = os.path.join(out_path, "busco_config.ini")
    with open(configfile, "w") as outF:
        config.write(outF)
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

# note: --force allows snakemake to handle rewriting files as necessary
# without needing to specify *all* busco outputs as snakemake outputs
shell(
    "busco --in {snakemake.input} --out {out_name} --force "
    " --cpu {snakemake.threads} --mode {mode} {lineage_cmd} "
    " {config_cmd} {extra} {log}"
)

# if not using a config to specific output location, move output.
# caveat: does not move output if user config incorrectly specifies out_path!
if out_path and not config_cmd:
    shell("cp -r {out_name} {out_path}")
    shell("rm -rf {out_name}")
