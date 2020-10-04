[![Snakemake](https://img.shields.io/badge/snakemake-≥5.18.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# as_analysis
A Snakemake pipeline for detecting allelic imbalance from DNA and RNA seq reads
![Pipeline Skeleton](https://drive.google.com/uc?export=view&id=1xefUeBPLLKfFfn9vu_IY2PHUBfC8XmST)

# download
Execute the following command.
```
git clone https://github.com/aryarm/as_analysis.git
```

# setup
The pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io). We recommend installing version 5.24.0:
```
conda create -n snakemake -c bioconda -c conda-forge --no-channel-priority 'snakemake==5.24.0'
```
We highly recommend you install [Snakemake via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) like this so that you can use the `--use-conda` flag when calling `snakemake` to let it [automatically handle all dependencies](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) of the pipeline. Otherwise, you must manually install the dependencies listed in the [env file](envs/default.yaml).

Note that our pipeline uses [WASP](https://github.com/bmvdgeijn/WASP) (v3.x), which is not available for download from conda. If you do not provide a local download of WASP, the pipeline will try to automatically install it.

# execution
1. Activate snakemake via `conda`:
    ```
    conda activate snakemake
    ```
2. Execute the pipeline

    Locally:
    ```
    ./run &
    ```
    __or__ on an SGE cluster:
    ```
    ./run --sge-cluster &
    ```
Log files describing the output of the pipeline will be created within the output directory. The `log` file contains a basic description of the progress of each rule, while the `qlog` file is more detailed.

### Executing the pipeline on your own data
You must modify [the config.yaml file](config.yaml) to specify paths to your data before you perform step 2 above. For more information about what is required in the config file, see the [READMEs for each portion of the pipeline](Snakefiles/README.md).

### Executing each portion of the pipeline separately
The entire pipeline is made up of three different sections. We provide a single [Snakefile](Snakefile) to execute all of them at once, but you can also execute each of these sections on their own. For each section that you'd like to run separately, you must fill out a new config file. You can find more information about these individual portions of the pipeline and how to execute them in the [Snakefiles directory](Snakefiles).

### Executing the pipeline on Google Cloud
```
./run-gcp &
```
See the [Google Cloud README](README.gcp.md) for full instructions.

### If this is your first time using Snakemake
We recommend that you run `snakemake --help` to learn about Snakemake's options. For example, to check that the pipeline will be executed correctly before you run it, you can call Snakemake with the `-n -p -r` flags. This is also a good way to familiarize yourself with the steps of the pipeline and their inputs and outputs (the latter of which are inputs to the first rule in each workflow -- ie the `all` rule).

Note that Snakemake will not recreate output that it has already generated, unless you request it. If a job fails or is interrupted, subsequent executions of Snakemake will just pick up where it left off. This can also apply to files that *you* create and provide in place of the files it would have generated.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# files
### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis. It uses DNA and RNA FASTQ files to generate a summary of allelic imbalance for each gene.

### config.yaml
Defines options and input for the Snakemake pipeline.

### run
An example bash script for executing the pipeline using `snakemake` and `conda`. Any arguments to this script are passed directly to `snakemake`.
