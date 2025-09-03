# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1)
This is the official repository for our paper [Genomic heterogeneity inflates the performance of variant pathogenicity predictions](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1).

It provides a **genome-wide, variant-type-stratified benchmark dataset** (>250,000 ClinVar variants) and the code to evaluate state-of-the-art DNA-based and protein-based models for variant pathogenicity prediction.

## Contents

- [Dataset](#dataset)
- [Notebook](#notebook)
- [Citation](#citation)

## Dataset
  
The dataset is available in the [`data/`](data/) directory.

- **Source**: ClinVar (GRCh38, single-nucleotide variants with expert-reviewed labels)  
- **Scale**: >250,000 variants  
- **Labels**: Pathogenic (1) vs. Benign (0) 
- **Variant categories**:  
  - Coding: missense, synonymous, start-loss, stop-gain, stop-loss  
  - Noncoding: canonical splice, non-splice intron, 5′ UTR, 3′ UTR, RNA gene

## Notebook

We provide one-click Jupyter notebook examples for each evaluated model.  

- **DNA-based models**:  
  AlphaGenome, DNABERT2, Evo2, GPN-MSA, Nucleotide Transformer (NT), PhyloGPN, PhyloP  
  Notebooks are available in the [`DNA-based Models/`](DNA-based%20Models/) directory.  

- **Protein-based models**:  
  ESM family models, AlphaMissense, PrimateAI-3D  
  Notebooks are available in the [`protein_models/`](protein_models/) directory.

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

