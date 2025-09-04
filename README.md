# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1)
This is the official repository for our paper [Genomic heterogeneity inflates the performance of variant pathogenicity predictions](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1).

It provides a **genome-wide, variant-type-stratified benchmark dataset** (>250,000 ClinVar variants) and the code to evaluate state-of-the-art DNA-based and protein-based models for variant pathogenicity prediction.

## Contents

- [Notebooks](#notebooks)
- [Results](#results)
- [Citation](#citation)

## Notebooks

We provide one-click Jupyter notebook examples for each evaluated model, benchmark creation, and results visualization.  

- **DNA-based models**:  
  AlphaGenome, DNABERT2, Evo2, GPN-MSA, Nucleotide Transformer (NT), PhyloGPN, PhyloP  
  → Notebooks are available in the [`DNA-based Models/`](DNA-based%20Models/) directory.  

- **Protein-based models**:  
  ESM family models, AlphaMissense, PrimateAI-3D  
  → Notebooks are available in the [`protein_models/`](protein_models/) directory.  

- **Benchmark creation**:  
  → See [`VEP_ClinVar_Benchmarking_RefSeq.ipynb`](VEP_ClinVar_Benchmarking_RefSeq.ipynb).  

- **Visualization**:  
  → See [`VEP_AUROC_figure.ipynb`](VEP_AUROC_figure.ipynb).


## Results

<p align="center">
  <img src="/figure1.png" alt="AUROC Results by Variant Type" width="900">
</p>

**Figure 1. Pathogenicity prediction performance of frontier sequence-based models across variant types.**  
Evaluation and comparison of DNA and protein sequence AI models for their capacity to distinguish between pathogenic and benign variants across variant types, measured by the area under the receiver operating characteristic curve (AUROC).  
- **%P** indicates the proportion of pathogenic variants in each group.  
- Some groups are defined by multiple annotated effects (e.g., both missense and 3′ UTR, with respect to different transcripts).  
- **DNA models** are shown as solid bars, **protein models** as dashed bars.  

*Note: The evaluation of PrimateAI-3D on stop-gain variants includes only 19,795 variants.*

## Citation

If you find this benchmark useful for your research, please cite our paper:

```bibtex
@article{YourName2025,
  author    = {Your Name and Collaborators},
  title     = {Genomic heterogeneity inflates the performance of variant pathogenicity predictions},
  journal   = {bioRxiv},
  year      = {2025},
  doi       = {10.1101/2025.XX.XXXXX},
  url       = {https://www.biorxiv.org/content/10.1101/2025.XX.XXXXXv1},
  eprint    = {https://www.biorxiv.org/content/10.1101/2025.XX.XXXXXv1.full.pdf}
}

