from dammit.meta import __path__

configfile: '../../databases.json'


rule download_databases:
    input: 'Pfam-A.hmm', 'Rfam.cm', 'OrthoDB.fasta', 'uniprot_sprot.fasta'


rule download_and_gunzip:
    output: '{database}.{file_type}'
    params:
        url = lambda wildcards: config[wildcards.database]['url'],
        md5 = lambda wildcards: config[wildcards.database].get('md5', False),
        metalink = lambda wildcards: config[wildcards.database].get('metalink', False)
    wrapper:
        f'file://{__path__}/wrappers/download'
        

