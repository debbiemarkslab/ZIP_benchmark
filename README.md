# ZIP_benchmark

Submitted to NeurIPS 2023 [Datasets and Benchmarks Track](https://neurips.cc/Conferences/2023/CallForDatasetsBenchmarks).

## Setup

1. Install the conda environment `indels_benchmark` as follows:

```python
conda env create --file indels_benchmark_conda.yml
```

The dependencies will take up about 4 GB of space if not already installed.

2. Download the required preprocessed datasets:

3. Install the baselines:
- [Progen2](https://github.com/salesforce/progen/tree/main/progen2)
- [Tranception](https://github.com/OATML-Markslab/Tranception) is included with minor changes as a subdirectory in this repo. Model checkpoints are available at TODO(Lood).

### Provean
To run Provean, the Provean software must be downloaded and installed from https://www.jcvi.org/research/provean#downloads. Note
that this also requires installing NCBI BLAST 2.4.0, CD-HIT 3.1.2 or greater, and a version 4 NCBI non-redundant protein database. In our experiments,
we use the non-redundant database at https://ftp.ncbi.nlm.nih.gov/blast/db/v4/. Download all files labeled "nr_*" to a folder, and then extract them.

### HMM
Follow the installation installation instructions in this fork of HMMER3: https://github.com/aaronkollasch/hmmer.
TODO(Aaron) 


### Optional:

#### VEP
If reproducing from raw files, the DDD-ASD dataset needs to be annotated with Ensembl-VEP. 
We used version 109, the latest version at time of publication (June 2023), available at http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html.

Note that the installer will download a cache, and the latest version of VEP caches and FASTA files will require a lot of storage (around 20GB). To avoid having to reannotate all the variants, download the preprocessed data.



## Reproducibility, data

The full processed benchmark dataset (except for the DDD dataset, which for legal reasons needs to be accessed and reconstructed separately) can be downloaded from (TODO Marks hosted dataset files) (X GB).

More information about the data sources and processing is provided in the Supplementary Materials of the paper (TODO link).

The `preprocessing` directory contains all code to reproduce the dataset from raw files. 

These raw files are available from (TODO Marks hosted raw files) and can be downloaded using the script in `preprocessing/download_raw_files.sh`.

Note: These raw files take up about X GB of space.


## Running baseline models
Scripts for running all of the baselines (Tranception, Progen-2, Provean, HMM) on each of the Datasets (ClinVar, DDD, DMS) is provided in the `baseline_models/bash_scripts/run_scripts/` directory.

After installing the conda environment, download the model checkpoints for ProGen2 and Tranception from here:

```TODO checkpoints```

## Reproducing plots
Code for reproducing plots is provided as Jupyter notebooks in the `figures` directory. (TODO upload these)

## Adding a baseline
To add a new baseline, add a a python script to the `baseline_models/scoring/` directory, and add a line running the model to each of the datasets' respective scripts in the `baseline_models/bash_scripts/run_scripts/` directory.

The input to the model will be a mapping file in CSV format, which contains a row per file of variants to be scored, along with other information such as the path to the alignment for each set of variants. Example file: `/processed_data/ClinVar/2023-05-05-tranception_mapping_by_refseq_id.csv`.

Each set of variants (e.g. an assay, or gene) is then containing one row per variant, with the mutant ID (e.g. the [HGVS annotation](http://varnomen.hgvs.org/recommendations/protein/) for clinical data), the wild-type protein sequence, and the variant protein sequence, with these column names provided in the mapping file. 

As an example, see `baseline_models/scoring/score_tranception.py`.

We'd love new baselines for comparison - please add a pull request on GitHub.

## Evaluation

The `stats/` directory contains scripts to evaluate models predictions - run the `stats/run_stats.sh` script to evaluate all of the baselines.

The clinical datasets (ClinVar-gnomAD, and DDD-ASD) use ROC-AUC and AUPRC binary metrics to evaluate the model's predictions distinguishing pathogenic from benign variants. Critically, these metrics are across all genes (not computed within genes), as there is such a low number of variants per gene.

The DMS assays have continuous outputs, so the Spearman rank correlation between the model's scores and the assay values are reported.

## Benchmark Datasets 
### gnomAD and ClinVar Pathogenicity Prediction 
### Deep Mutational Scans 
### DDD/ASD Developmental Disorder Prediction

## Multiple Sequence Alignments
Tranception with Retrieval uses a Multiple Sequence Alignment (MSA) as part of its input. We generate MSAs for each 
of our new datasets and make them available for download here (add link later). The details of how the alignments were generated are described below:

### gnomAD and ClinVar 

### DMS experiments 

### DDD and ASD

## Contributing

We welcome new baselines (see above) as well as general improvement pull requests.

## Citation

Please cite this benchmark if you use it in your work:
TODO
