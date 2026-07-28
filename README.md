# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2)

This repository contains the scoring and figure-generation notebooks for our paper, [Genomic heterogeneity inflates the performance of variant pathogenicity predictions](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2).

All notebooks are stand-alone and use repository-relative paths. Model-scoring notebooks can be opened from the repository root or directly from their model subdirectory; figure notebooks should be run from the repository root. Comments and descriptions are in English.

The repository includes [`data/sample_data.csv.gz`](data/sample_data.csv.gz), a representative 1,000-row sample containing the annotations needed by the model-scoring notebooks. It allows readers to test the scoring code without first downloading the complete benchmark.

Figure notebooks use the complete scored benchmark at `data/processed/clinvar_benchmark_updated.csv`. This file is intentionally excluded from GitHub because of its size. The notebooks retain the outputs from the full-data analysis reported in the manuscript.

## Repository contents

### DNA-based model scoring

The [`DNA-based Models/`](DNA-based%20Models/) directory contains notebooks for:

- AlphaGenome, including the all-output analysis used for Supplementary Figure S5
- DNABERT-2, pinned to model revision `7bce263b15377fc15361f52cfab88f8b586abda0`
- Evo 2
- GPN-MSA
- GPN-Star
- Nucleotide Transformer v3
- PhyloGPN, pinned to the April 2026 model revision `3556db4c469e67d25f0f7a0a6653b48be3eebf51`
- PhyloP
- Rule-based baseline

### Protein-based model scoring

The [`protein_models/`](protein_models/) directory contains notebooks for:

- AlphaMissense
- ESM1b, ESM1v, and ESM2
- PrimateAI-3D (the [official hg38 score table](https://huggingface.co/datasets/illumina-ai/PrimateAI-3D-scores) is gated by its academic license)
- VESM++

### Benchmark construction

[`VEP_ClinVar_Benchmarking_RefSeq.ipynb`](VEP_ClinVar_Benchmarking_RefSeq.ipynb) builds the annotation-only ClinVar benchmark.

Its fixed public inputs are:

- [ClinVar variant summary, February 2026 archive](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-02.txt.gz)
- [MANE Select GRCh38 release 1.5](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.refseq_genomic.gff.gz)
- [UCSC hg38 reference genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- [UCSC NCBI RefSeq transcript table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz)

### Figure generation

| Manuscript figure | Notebook |
|---|---|
| Figure 1: variant-type AUROC | [`VEP_AUROC_figure.ipynb`](VEP_AUROC_figure.ipynb) |
| Figure 2: model robustness | [`VEP_model_robustness_figure.ipynb`](VEP_model_robustness_figure.ipynb) |
| Figure 3: score distributions | [`VEP_score_distribution_figure.ipynb`](VEP_score_distribution_figure.ipynb) |
| Supplementary Figure S2: coordinate-balanced AUROC | [`VEP_coordinate_balance_figure.ipynb`](VEP_coordinate_balance_figure.ipynb) |
| Supplementary Figure S3: ClinVar review-star analysis | [`VEP_ClinVar_star_figure.ipynb`](VEP_ClinVar_star_figure.ipynb) |
| Supplementary Figure S4: AUPRC | [`VEP_AUPRC_figure.ipynb`](VEP_AUPRC_figure.ipynb) |
| Supplementary Figure S5: AlphaGenome splice outputs | [`VEP_AlphaGenome_splice_figure.ipynb`](VEP_AlphaGenome_splice_figure.ipynb) |

## Input and output

Model-scoring notebooks use the included sample by default:

```text
data/sample_data.csv.gz
```

Figure notebooks use the complete scored benchmark:

```text
data/processed/clinvar_benchmark_updated.csv
```

Model-scoring notebooks write new scores to:

```text
data/model_scores/
```

Figure notebooks write generated tables and figures to:

```text
results/tables/
results/figures/
```

The full benchmark and generated model-score files are intentionally excluded from GitHub. The stable Supplementary Data download URL will be added here when it becomes available.

## Running the notebooks

For figure generation, install the common analysis packages:

```bash
python -m pip install jupyter pandas numpy scipy scikit-learn matplotlib seaborn
jupyter lab
```

Model-scoring dependencies differ by model and are stated inside each model notebook. Several scoring notebooks require a GPU, gated model access, an API key, or a manually licensed dataset.

## Results

<p align="center">
  <img src="/Figure1.svg" alt="AUROC results by variant type" width="900">
</p>

## Citation

```bibtex
@article{genomic2025biorxiv,
  author  = {Baiyu Lu and Xueshen Liu and Po-Yu Lin and Nadav Brandes},
  title   = {Genomic heterogeneity inflates the performance of variant pathogenicity predictions},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.09.05.674459},
  url     = {https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2}
}
```

## License

Code in this repository is available under the [MIT License](LICENSE). External datasets and model weights retain their original licenses.
