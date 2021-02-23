[![Build Status](https://travis-ci.org/mpozuelo/MGI_demux.svg?branch=master)](https://travis-ci.org/mpozuelo/MGI_demux)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/rnaseq.svg)](https://hub.docker.com/r/nfcore/rnaseq/)

### Introduction

**mpozuelo/index_BC** is a bioinformatics analysis pipeline used for getting occurrences at positions for i7 (index), i5(index2) and bc (Rd1) at CIMA cluster.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Start running your own analysis!

Run in a linux/MACOs machine and there is no need to use singularity nor docker

```bash
nextflow run mpozuelo/index_BC -profile <docker/singularity/conda> --input 'samplesheet.txt'
```

def run = lane.run
def lane = lane.lane
def bcsize = lane.bc
def indexsize = lane.index
def index2size = lane.index2
def libsize = lane.size
The samplesheet used must have at least 6 columns comma separated WITH header: run, lane, index, index2, bc, size

size: 100 (PE100), 200 (PE200)
index, index2 and bc: number indicating the size of the sequence for each of these barcodes

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

### Documentation

The mpozuelo/bn_index pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Running the pipeline](docs/usage.md)

### Credits

These scripts were originally written for getting occurrences at positions for i7 (index), i5(index2) and bc (Rd1) at CIMA cluster by Marta Pozuelo ([@mpozuelo](https://github.com/mpozuelo))
