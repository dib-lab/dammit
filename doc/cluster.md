# Distributing dammit jobs across a cluster

`dammit` can run on a single compute instance, or can submit each individual job to a job scheduler, 
if you provide the right submission information for your cluster. Job submission is handled
via snakemake, so please see the [snakemake cluster documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
for the most up-to-date version of these instructions.

## Using A Snakemake Profile for Job Submission

### Set up a snakemake profile for your cluster

We recommend using a [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
to enable job submission for the job scheduler used by your cluster. You can start from the cookiecutter 
profiles [here](https://github.com/snakemake-profiles/doc) or write your own.

### Direct dammit to use the snakemake profile 

When you'd like dammit to submit jobs to a job scheduler, direct it to use your cluster profile by 
adding `--profile <profile-folder-name>` at or near the end of your dammit command (after all dammit-specific arguments).
Again, see the [snakemake profile documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for additional information.

