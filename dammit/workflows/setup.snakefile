import os
from dammit.meta import __path__


rule run_dammit_rename_transcriptome:
    input:
        config["input_transcriptome"],
    output:
        fasta=os.path.join(results_dir, "{transcriptome}.fasta"),
        names=os.path.join(results_dir, "{transcriptome}.namemap.csv")
    log:
        os.path.join(logs_dir, "rename-transcriptome.log")
    params:
        basename = config.get("basename", "Txome"),
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'
