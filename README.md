# P3 Model Iteration

This repo contains code for a simple Snakemake-based pipeline focussed on testing out
various modeling approaches and feature and response data variants.

## Usage

To setup the pipeline, clone the repo and from inside the repo directory, run:

```
conda create -n p3-models --file requirements.txt
conda activate p3-models
```

This will create a new conda environment containing all of the dependencies needed by
the pipline.

Next, edit the configuration file (`config.yml`) to indicate the desired settings and
locations of expected input files.

To verify that the pipeline is configured properly and is able to locate the neccessary
files, you can perform a dry run using:

```
snakemake -n
```

To run the pipeline, simply calling Snakemake, specifying the desired number of threads,
e.g. 

```
snakemake -j 8
```
