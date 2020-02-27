__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

protein_cmd = ""
protein = snakemake.params.get("protein_input", False)

if protein == "prot":
    protein = True
elif protein == "nucl":
    protein = False

if protein:
    protein_cmd = " -p "

out_base = str(snakemake.output)[:-4]

shell("lastdb {extra} {protein_cmd} -v -P {snakemake.threads} {out_base} {snakemake.input} {log}")
