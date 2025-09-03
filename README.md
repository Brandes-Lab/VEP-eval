# Genomic heterogeneity inflates the performance of variant pathogenicity predictions

[![bioRxiv](https://img.shields.io/badge/bioRxiv-Available-red)](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1)
This is the official repository for our paper [Genomic heterogeneity inflates the performance of variant pathogenicity predictions](https://www.biorxiv.org/content/XXXX.XX.XXXXXXv1).

It provides a **genome-wide, variant-type-stratified benchmark dataset** (>250,000 ClinVar variants) and the code to evaluate state-of-the-art **DNA- and protein-based models** for variant pathogenicity prediction.

## Contents

- [Dataset](#dataset)
- [Notebook](#notebook)
- [Results](#results)
- [Citation](#citation)

## Dataset
  
The dataset is available in the [`data/`](data/) directory.

- **Source**: ClinVar (GRCh38, single-nucleotide variants with expert-reviewed labels)  
- **Scale**: >250,000 variants  
- **Labels**: Pathogenic vs. Benign  
- **Variant categories**:  
  - Coding: missense, synonymous, start-loss, stop-gain, stop-loss  
  - Noncoding: canonical splice, non-splice intron, 5′ UTR, 3′ UTR, RNA gene


This repository contains all the codes to replicate the benchmark generation, the scores of protein sequence and DNA sequence mutations. <br>
We categorize the model into 2 directories: <br>
DNA models: AlphaGenome, DNABERT2, Evo2, GPNMSA, NT, PhyloGPN, PhyloP. <br>
Protein models: ESM Models, AlphaMissense, PrimateAI3D. <br>
