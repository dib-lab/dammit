"""Snakemake wrapper for Transdecoder Predict"""

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output_dir = os.path.dirname(str(snakemake.input.longorfs))

addl_outputs = ""
pfam = snakemake.input.get("pfam_hits", "")
if pfam:
    addl_outputs += " --retain_pfam_hits " + pfam

blast = snakemake.input.get("blastp_hits", "")
if blast:
    addl_outputs += " --retain_blastp_hits " + blast

input_fasta = str(snakemake.input.fasta)
if input_fasta.endswith("gz"):
    input_fa = input_fasta.rsplit(".gz")[0]
    shell("gunzip -c {input_fasta} > {input_fa}")
else:
    input_fa = input_fasta

shell("TransDecoder.Predict --output_dir {output_dir} -t {input_fa} {addl_outputs} {extra} {log}")

outputs = snakemake.output
for outp in outputs:
    old_output = os.path.basename(outp)
    shell("mv {old_output} {outp}")

shell("mv pipeliner* {output_dir}")
