__author__ = "Camille Scott and N. Tessa Pierce"
__copyright__ = "Copyright 2019, Camille Scott"
__email__ = "cswel@ucdavis.edu"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.params.get('metalink') and snakemake.params.get('md5'):
    raise ValueError('Do not use metalink and md5 simultaneously.')

cmd = ['curl', '-L']

cmd.append('{snakemake.params.url}')

if snakemake.params.get('metalink'):
    cmd.extend(['--metalink', '{snakemake.params.metalink}'])

donefile = None
output =  str(snakemake.output)
if output.endswith(".done"):
    donefile = output
    output = output.rsplit(".done")[0]

fileformat = snakemake.params.get("fileformat")
if fileformat == "gz":
    cmd.extend(['|', 'gunzip', '-c', '>', '{output}'])
elif fileformat == "tar.gz":
    cmd.extend(['|', 'tar', '-xzf', '>', '{output}'])
else:
    raise ValueError('Valid filetypes are "gz" and "tar.gz"')

if snakemake.params.get('md5'):
    cmd.append('&& python -c "assert \'`md5sum {output} | '
                'awk \'{{print $1}}\'`\' == \'{snakemake.params.md5}\', '
                '\'MD5sum does not match\'"')

cmd.append("{log}")
shell(' '.join(cmd))

if donefile:
    shell("touch {donefile}")
