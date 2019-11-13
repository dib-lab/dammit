__author__ = "Camille Scott"
__copyright__ = "Copyright 2019, Camille Scott"
__email__ = "cswel@ucdavis.edu"
__license__ = "MIT"

from snakemake.shell import shell

if snakemake.params.get('metalink') and snakemake.params.get('md5'):
    raise ValueError('Do not use metalink and md5 simultaneously.')

cmd = ['curl', '-L']

if snakemake.params.get('metalink'):
    cmd.extend(['--metalink', '{snakemake.params.metalink}'])

cmd.append('{snakemake.params.url}')
cmd.extend(['|', 'gunzip', '-c', '>', '{snakemake.output}'])

if snakemake.params.get('md5'):
    cmd.append('&& python -c "assert \'`md5sum {snakemake.output} | '
                'awk \'{{print $1}}\'`\' == \'{snakemake.params.md5}\', '
                '\'MD5sum does not match\'"')

print(' '.join(cmd))
shell(' '.join(cmd))
