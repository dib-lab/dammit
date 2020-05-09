""" Snakemake wrapper for shmlast """

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2020, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

evalue = snakemake.params.get("evalue", "0.000001")

search_type = snakemake.params.get("search_type", "crbl")
assert search_type in ["rbl", "crbl"], "search_type must be either 'rbl' or 'crbl'"


shell(
    "shmlast {search_type} -q {snakemake.input.query} -d {snakemake.input.database} --n_threads {snakemake.threads} --evalue-cutoff {evalue} -o {snakemake.output} {log}"
)
