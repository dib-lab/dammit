import os
from dammit.meta import __path__


rule dammit_rename_transcriptome:
    input:
        config["input_transcriptome"],
    output:
        fasta = os.path.join(results_dir, "{transcriptome}.fasta"),
        names = os.path.join(results_dir, "{transcriptome}.namemap.csv")
    log:
        os.path.join(logs_dir, "{transcriptome}.rename.log")
    params:
        basename = config.get("basename", "Txome"),
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'
