"""Snakemake wrapper for hmmscan"""

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

profile = snakemake.input.get("profile")

profile = profile.rsplit(".h3", 1)[0]
assert profile.endswith(".hmm"), 'your profile file should end with ".hmm" '

# Direct the main human-readable output to a file <f> instead of the default stdout.
out_cmd = ""
outfile = snakemake.output.get("outfile", "")
if outfile:
    out_cmd += " -o {} ".format(outfile)

# save parseable table of per-domain hits to file <f>
domtblout = snakemake.output.get("domtblout", "")
if domtblout:
    out_cmd += " --domtblout {} ".format(domtblout)
else:
    assert False, "must specificy domtblout in outputs"

# default params: enable evalue threshold. If bitscore thresh is provided, use that instead (both not allowed)
# report models >= this score threshold in output
evalue_threshold = snakemake.params.get("evalue_threshold", 0.00001)
score_threshold = snakemake.params.get("score_threshold", "")

if score_threshold:
    thresh_cmd = " -T {} ".format(float(score_threshold))
else:
    thresh_cmd = " -E {} ".format(float(evalue_threshold))

# all other params should be entered in "extra" param
extra = snakemake.params.get("extra", "")

if snakemake.threads == 1:
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    shell(
        "hmmscan {out_cmd} {thresh_cmd} --cpu {snakemake.threads}"
        " {extra} {profile} {snakemake.input.fasta} {log}"
    )
else:
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)
    shell(
        "ope parallel -j {snakemake.threads} {snakemake.input.fasta} "
        "hmmscan --cpu 1 -o /dev/null --domtblout /dev/stdout {thresh_cmd} {extra} {profile} /dev/stdin > {domtblout} {log}"
    )
