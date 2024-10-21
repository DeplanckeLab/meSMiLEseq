This repository contains code to reproduce analyses presented in the paper 'Identification of methylation-sensitive human transcription factors using meSMiLE-seq'. 

> [!NOTE]
> **Prerequisites**
> 
> A prerequisite to running the code provided in this repo is obtaining the rights to use the **ProBound Suite** (Rube et al., 'Prediction of protein-ligand binding affinity from sequencing data with interpretable machine learning', doi:10.1038/s41587-022-01307-0) and a **material transfer agreement**. The *pyProBound operator package* must be installed to use the ProBound Suite via Python. The package source code can be found at **probound_operator**.

## Structure

Each folder consists of subfolders containing exemplary data, scripts, and an output folder. The analysis can either be executed as a full block by calling 'run_pipeline.sh' or it can be split into separate parts.

**Contents**
1. **SMiLEseq**
- analysis of classical SMiLE-seq experiments by using a Fisher's exact test to prefilter raw sequencing reads before *de novo* motif discovery via ProBound

2. **meSMiLEseq**
- analysis of methylation-sensitive SMiLE-seq data to create methylation-aware binding models and *k*-mer scatterplots

3. **WGBS_analysis**
- exemplary analysis pipeline using PRDM13 to show CG methylation patterns at individual motif occurrences in cells by intersecting ChIP-seq and WGBS data

4. **Motifs**
- folder contains all TF binding motifs showcased in **Figure 1** as 'position-specific affinity matrices' and DNA logos

5. **probound_operator**
- folder contatins source code of a package to operate ProBound from Python, provided ProBound is installed and its absolute path is stored in environment variable PROBOUND_JAR_FULL_PATH. You can create the package using the following code:
```
export PROBOUND_JAR_FULL_PATH="path/to/probound/jar"
cd probound_operator
python3 -m build
```
