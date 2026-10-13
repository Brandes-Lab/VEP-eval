# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![Nature Communications](https://img.shields.io/badge/Nature%20Communications-2026-blue)](https://doi.org/10.1038/s41467-026-78623-5)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2)
[![Data](https://img.shields.io/badge/Hugging%20Face-dataset-yellow)](https://huggingface.co/datasets/brandeslab/VEP-eval)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22716331-blue)](https://doi.org/10.5281/zenodo.22716331)

We provide the scoring, benchmark-construction, and figure-generation notebooks for our paper, *Genomic heterogeneity inflates the performance of variant pathogenicity predictions* (Nature Communications, 2026; [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.05.674459v2)).

We designed the notebooks to use repository-relative paths and English comments throughout. Model-scoring notebooks can be opened from the repository root or directly from their model subdirectory; benchmark and figure notebooks should be run from the repository root.

We include [`data/sample_data.csv.gz`](data/sample_data.csv.gz), a representative 1,000-row sample containing the annotations required by the model-scoring notebooks. This sample allows the scoring code to be tested without first downloading the complete benchmark.

The figure notebooks use the complete scored benchmarks, which we host on [Hugging Face](https://huggingface.co/datasets/brandeslab/VEP-eval) rather than on GitHub because of their size. `python scripts/download_data.py` downloads them and places them where the notebooks expect them (see [Reproducing the figures](#reproducing-the-figures)). The executed figure notebooks retain the outputs from the complete-data analysis.

## Repository contents

### DNA-based model scoring

In [`DNA-based Models/`](DNA-based%20Models/), we provide scoring notebooks for:

- AlphaGenome, including the all-output analysis used for Supplementary Figure 6
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

We use [`VEP_ClinVar_Benchmarking_RefSeq.ipynb`](VEP_ClinVar_Benchmarking_RefSeq.ipynb) to build the annotation-only ClinVar benchmark, following the filtering workflow summarized in Supplementary Figure 1.

Its fixed public inputs are:

- [ClinVar variant summary, February 2026 archive](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-02.txt.gz)
- [MANE Select GRCh38 release 1.5](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.refseq_genomic.gff.gz)
- [UCSC hg38 reference genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- [UCSC NCBI RefSeq transcript table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz)

We use [`VEP_COSMIC_Benchmarking.ipynb`](VEP_COSMIC_Benchmarking.ipynb) to build the somatic-variant benchmark from COSMIC Cancer Mutation Census v103. The notebook gives the exact local filenames and download command for MANE Select. CMC v103 must be downloaded after login from the [official COSMIC Cancer Mutation Census page](https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census); we do not redistribute the licensed COSMIC source file.

### Figure generation

| Manuscript figure | Notebook | Input | Approx. runtime |
|---|---|---|---|
| Figure 1: variant-type AUROC | [`VEP_AUROC_figure.ipynb`](VEP_AUROC_figure.ipynb) | ClinVar benchmark | 35 min |
| Figure 2: splice and 5′ UTR score distributions | [`VEP_score_distribution_figure.ipynb`](VEP_score_distribution_figure.ipynb) | ClinVar benchmark | < 1 min |
| Figure 3: model robustness across variant types | [`VEP_model_robustness_figure.ipynb`](VEP_model_robustness_figure.ipynb) | ClinVar benchmark | < 1 min |
| Supplementary Figure 1: ClinVar benchmark construction flowchart | [`VEP_benchmark_pipeline_figure.ipynb`](VEP_benchmark_pipeline_figure.ipynb) | none | < 1 min |
| Supplementary Figure 2: coordinate-balanced AUROC | [`VEP_coordinate_balance_figure.ipynb`](VEP_coordinate_balance_figure.ipynb) | ClinVar benchmark | 13 min |
| Supplementary Figure 3: ClinVar review-star analysis | [`VEP_ClinVar_star_figure.ipynb`](VEP_ClinVar_star_figure.ipynb) | ClinVar benchmark | 22 min |
| Supplementary Figure 4: minority-class AUPRC | [`VEP_AUPRC_figure.ipynb`](VEP_AUPRC_figure.ipynb) | ClinVar benchmark | 23 min |
| Supplementary Figure 5: COSMIC somatic-variant AUROC | [`VEP_COSMIC_AUROC_figure.ipynb`](VEP_COSMIC_AUROC_figure.ipynb) | COSMIC benchmark | 7 min |
| Supplementary Figure 6: AlphaGenome splice outputs | [`VEP_AlphaGenome_splice_figure.ipynb`](VEP_AlphaGenome_splice_figure.ipynb) | ClinVar benchmark | < 1 min |

Runtimes are for one notebook on a laptop CPU; the 1,000-resample stratified bootstraps dominate. Each notebook is independent and can be run on its own.

## Input and output

Model-scoring notebooks use the included sample by default:

```text
data/sample_data.csv.gz
```

The figure notebooks read two complete benchmarks, written there by `scripts/download_data.py`:

```text
data/processed/clinvar_benchmark_updated.csv   # ClinVar benchmark, 242,132 variants
cosmic/cosmic_all_scores.csv                   # COSMIC benchmark, 55,538 variants (Supplementary Figure 5)
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

We exclude the complete benchmarks and generated model-score files from GitHub because of their size. The benchmarks, including every model score used in the paper, are available on [Hugging Face](https://huggingface.co/datasets/brandeslab/VEP-eval).

## Reproducing the figures

Run every command from the repository root.

1. Create the environment. We generated the published figures with Python 3.11 and the package versions pinned in [`requirements.txt`](requirements.txt).

   ```bash
   python -m venv .venv && source .venv/bin/activate
   python -m pip install -r requirements.txt
   ```

2. Make sure the Arial font is installed. All figure text is set in Arial; without it matplotlib substitutes another font and the layout changes. Arial ships with macOS and Windows; on Debian or Ubuntu, install `ttf-mscorefonts-installer` and then run `rm -rf ~/.cache/matplotlib`.

3. Download the benchmarks. The script fetches them from a pinned Hugging Face revision and verifies their SHA-256 checksums.

   ```bash
   python scripts/download_data.py
   ```

4. Execute the figure notebooks.

   ```bash
   for nb in VEP_*_figure.ipynb; do
     jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 "$nb"
   done
   ```

   Figures are written to `results/figures/` as PDF, SVG and PNG, and the values behind them to `results/tables/`.

With this environment, the notebooks reproduce the published Figures 1–3 and Supplementary Figures 2–6 exactly: we checked every plotted value against the regenerated tables and compared the output artwork with the published figures. Supplementary Figure 1 is a schematic of the benchmark filtering steps; its counts are entered directly in `VEP_benchmark_pipeline_figure.ipynb`.

Model-scoring dependencies differ by model and are stated inside each notebook. Several models require a GPU, gated model access, an API key, or a separately licensed dataset.

## Results

<p align="center">
  <img src="/Figure1.svg" alt="AUROC results by variant type" width="900">
</p>

## Citation

```bibtex
@article{lu2026genomic,
  author  = {Lu, Baiyu and Liu, Xueshen and Lin, Po-Yu and Brandes, Nadav},
  title   = {Genomic heterogeneity inflates the performance of variant pathogenicity predictions},
  journal = {Nature Communications},
  year    = {2026},
  doi     = {10.1038/s41467-026-78623-5}
}
```

To cite this code, please use the Zenodo record [10.5281/zenodo.22716331](https://doi.org/10.5281/zenodo.22716331), which always resolves to the latest release.

## License

Code in this repository is available under the [MIT License](LICENSE). External datasets and model weights retain their original licenses.
