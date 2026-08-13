# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2)

We provide the scoring, benchmark-construction, and figure-generation notebooks for our paper, [Genomic heterogeneity inflates the performance of variant pathogenicity predictions](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2).

We designed the notebooks to use repository-relative paths and English comments throughout. Model-scoring notebooks can be opened from the repository root or directly from their model subdirectory; benchmark and figure notebooks should be run from the repository root.

We include [`data/sample_data.csv.gz`](data/sample_data.csv.gz), a representative 1,000-row sample containing the annotations required by the model-scoring notebooks. This sample allows the scoring code to be tested without first downloading the complete benchmark.

Our ClinVar figure notebooks use the complete scored benchmark at `data/processed/clinvar_benchmark_updated.csv`. We do not store this file on GitHub because of its size. The executed figure notebooks retain the outputs from the complete-data analysis.

## Repository contents

### DNA-based model scoring

In [`DNA-based Models/`](DNA-based%20Models/), we provide scoring notebooks for:

- AlphaGenome, including the all-output analysis used for Supplementary Figure S6
- DNABERT-2, pinned to model revision `7bce263b15377fc15361f52cfab88f8b586abda0`
- Evo 2
- GPN-MSA
- GPN-Star
- Nucleotide Transformer v3
- PhyloGPN, pinned to the April 2026 model revision `3556db4c469e67d25f0f7a0a6653b48be3eebf51`
- PhyloP
- Rule-based baseline

### Protein-based model scoring

In [`protein_models/`](protein_models/), we provide scoring notebooks for:

- AlphaMissense
- ESM1b, ESM1v, and ESM2
- PrimateAI-3D (the [official hg38 score table](https://huggingface.co/datasets/illumina-ai/PrimateAI-3D-scores) is gated by its academic license)
- VESM++

### Benchmark construction

We use [`VEP_ClinVar_Benchmarking_RefSeq.ipynb`](VEP_ClinVar_Benchmarking_RefSeq.ipynb) to build the annotation-only ClinVar benchmark and the workflow summarized in Supplementary Figure S1.

Its fixed public inputs are:

- [ClinVar variant summary, February 2026 archive](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-02.txt.gz)
- [MANE Select GRCh38 release 1.5](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.refseq_genomic.gff.gz)
- [UCSC hg38 reference genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- [UCSC NCBI RefSeq transcript table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz)

We use [`VEP_COSMIC_Benchmarking.ipynb`](VEP_COSMIC_Benchmarking.ipynb) to build the somatic-variant benchmark from COSMIC Cancer Mutation Census v103. The notebook gives the exact local filenames and download command for MANE Select. CMC v103 must be downloaded after login from the [official COSMIC Cancer Mutation Census page](https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census); we do not redistribute the licensed COSMIC source file.

### Figure generation

| Manuscript figure | Notebook |
|---|---|
| Figure 1: variant-type AUROC | [`VEP_AUROC_figure.ipynb`](VEP_AUROC_figure.ipynb) |
| Figure 2: splice and 5′ UTR score distributions | [`VEP_score_distribution_figure.ipynb`](VEP_score_distribution_figure.ipynb) |
| Figure 3: model robustness across variant types | [`VEP_model_robustness_figure.ipynb`](VEP_model_robustness_figure.ipynb) |
| Supplementary Figure S1: ClinVar benchmark workflow | [`VEP_ClinVar_Benchmarking_RefSeq.ipynb`](VEP_ClinVar_Benchmarking_RefSeq.ipynb) |
| Supplementary Figure S2: coordinate-balanced AUROC | [`VEP_coordinate_balance_figure.ipynb`](VEP_coordinate_balance_figure.ipynb) |
| Supplementary Figure S3: ClinVar review-star analysis | [`VEP_ClinVar_star_figure.ipynb`](VEP_ClinVar_star_figure.ipynb) |
| Supplementary Figure S4: minority-class AUPRC | [`VEP_AUPRC_figure.ipynb`](VEP_AUPRC_figure.ipynb) |
| Supplementary Figure S6: AlphaGenome splice outputs | [`VEP_AlphaGenome_splice_figure.ipynb`](VEP_AlphaGenome_splice_figure.ipynb) |

## Input and output

Model-scoring notebooks use the included sample by default:

```text
data/sample_data.csv.gz
```

Our ClinVar figure notebooks use the complete scored benchmark:

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

We exclude the complete benchmarks and generated model-score files from GitHub because of file size and, for COSMIC, source-data licensing. We will add the stable Supplementary Data URL here when it becomes available.

## Running the notebooks

To run the figure notebooks, install the common analysis packages:

```bash
python -m pip install jupyter pandas numpy scipy scikit-learn matplotlib seaborn
jupyter lab
```

Model-scoring dependencies differ by model and are stated inside each notebook. Several models require a GPU, gated model access, an API key, or a separately licensed dataset.

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
