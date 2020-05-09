import os
from dammit.meta import __path__


rule dammit_rename_transcriptome:
    input:
        config["input_transcriptome"],
    output:
        fasta=os.path.join(results_dir, "{transcriptome}.fasta"),
        names=os.path.join(results_dir, "{transcriptome}.namemap.csv")
    log:
        os.path.join(logs_dir, "{transcriptome}.rename.log")
    params:
        basename = config.get("basename", "Txome"),
    script: f'file://{__path__}/wrappers/dammit/rename-transcriptome.wrapper.py'

#rule back_to_original_names:
#    input:
#        fasta=os.path.join(results_dir, "{transcriptome}.dammit.fasta"),
#        namemap=os.path.join(results_dir, "{transcriptome}.dammit.namemap"),
#        gff3=os.path.join(results_dir, "{transcriptome}.dammit.gff3")
#    output:
#        fasta=os.path.join(results_dir, "{transcriptome}.orignames.fasta"),
#        gff3=os.path.join(results_dir, "{transcriptome}.orignames.fasta")
#    log:
#        os.path.join(logs_dir, "{transcriptome}.convert_to_originalnames.log")
#    benchmark:
#        os.path.join(logs_dir, "{transcriptome}.convert_to_originalnames.benchmark")
#    shell:
#        """
#        python dammit_gff3_to_trinity_names.py {input.fasta} --dammit_namemap {input.namemap}
#        --dammit_gff3 {input.gff3} --outFasta {output.fasta} --outGFF {output.gff}
#        """
