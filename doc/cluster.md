# Distributing dammit jobs across a cluster

dammit can run on a single compute instance, or can submit each individual job to a job scheduler, 
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


## Resource Modification for Cluster Submission


### Limit maximum memory utilization

Each step of dammit runs with pre-defined default resources. You can limit the maximum amount of memory utilized at at a single time by providing the `--resources`
snakemake argument near the end of your dammit command (after all dammit-specific arguments). To limit memory, pass in the maximum memory you'd like the workflow
to use (in megabytes) --e.g. `--resources mem_mb=20000` limits dammit to `20Gb` of memory.


### Modify tool-specific memory allocation (Advanced Use Only!!)

To modify the memory limitation for jobs submitted for each specific step, you can pass in a custom configuration file
that modifies the resources specified to each snakemake rule within dammit. 

_[todo: add details if/when we enable this]_

