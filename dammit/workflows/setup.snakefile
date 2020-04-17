import os
from dammit.meta import __path__


rule run_dammit_rename_transcriptome:
    input:
        config["input_transcriptome"]
    output:
        os.path.join(results_dir, "{transcriptome}.fasta")
    log:
        os.path.join(logs_dir, "rename-transcriptome.log")
    params:
        basename = "Txome"
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'
