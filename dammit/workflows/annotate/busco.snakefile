import os
from dammit import meta

rule run_busco:
    input: {assembly} 
    output:
        "txome_busco/full_table_txome_busco.tsv",
    log:
        "logs/transcriptome_busco.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage_path=os.path.join(DATA_DIR, "example-busco-lineage"),
        # optional parameters
        extra="" 
    conda: "environment.yml"
    script: "busco.wrapper.py"

#def parse_busco_summary:
