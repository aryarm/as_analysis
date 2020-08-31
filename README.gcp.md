# Executing our pipeline on Google Cloud
We offer preliminary cloud support for the WASP and counts pipelines. Follow these instructions to run those pipelines on GTEX data.

## Setup
This setup generally follows [Snakemake's Google Life Sciences Tutorial](https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html).
1. Setup the [config file](configs/config-WASP.yaml)
2. Add your [Google Application Credentials and enable the requisite APIs](https://snakemake.readthedocs.io/en/stable/executing/cloud.html#executing-a-snakemake-workflow-via-google-cloud-life-sciences)
3. Create a cloud storage bucket and upload the following from your [config file](configs/config-WASP.yaml)
	- the `samples` file
	- the `chrom_info` file
	- the `gene_info` file
	- the `ref_genome` file
4. Use your ERA Commons Account to get access to the GTEX data through a Terra-based Google Cloud storage bucket
5. Copy the following GTEX data from the Terra bucket to your own:
    - the VCF and its `.tbi` index (to the path in your [config file](configs/config-WASP.yaml))
    - the BAM samples and their `.bai` index (to the [`map1_sort`](/Snakefiles/README.WASP.md#output) folder within the [config file](configs/config-WASP.yaml)'s `output_dir`)
6. Run the pipeline!
```
./run-gcp &
```

## Caveats
Some of the steps within our pipeline will not run properly on GCP. We've refrained from changing our pipeline to accommodate these problems because they are mostly related to bugs within Snakemake (or other things that we expect to improve with time).

issue | affected rules | workaround
---|---|---
[directory() output](https://github.com/snakemake/snakemake/issues/576) | create_STAR_index (and all downstream steps) | perform this step manually on your cluster, then reupload to the storage bucket
[checkpoints](https://github.com/snakemake/snakemake/issues/574) | vcf_chroms, split_vcf_by_chr, vcf2h5 | download the VCF to your cluster, perform the affected steps manually, then reupload to the storage bucket
[WASP](https://github.com/aryam7/WASP/issues/16) | get_WASP, install_WASP | add WASP to your working directory (under git version control)

## Other challenges


## [Planned features](https://github.com/aryam7/as_analysis/issues?q=is%3Aopen+is%3Aissue+label%3Agcp)
