# tcrBLOSUM: an amino acid substitution matrix for sensitive alignment of distant epitope-specific TCRs
tcrBLOSUM is a specialized BLOSUM-style amino acid substitution matrix tailored for TCR sequences. This matrix captures the amino acid variations occurring within TCRs binding the same epitope. Employing tcrBLOSUM instead of BLOSUM62 in TCR clustering methods like [TCRdist](https://github.com/kmayerb/tcrdist3) enables epitope-annotation of distant TCRs, with diverse amino acid composition and physicochemical profiles.

![Alt text](./results/figures/graph_abstract.png?raw=true "graphical abstract")

[//]: # (<p align="center">)

[//]: # (  <img src="results/figures/graph_abstract.png" alt="graphical abstract" width="800" />)

[//]: # (</p>)

## Overview

- *data/*: data used to construct **tcrBLOSUM** and two descriptor-based amino acid substitution matrices (**PhysChemSim** and **TopoSim**)
  - *data/parsed/processed_for_tcr_blosum_calculations*: TCR data with defined blocks of epitope-specific TCRs

- *src/*: scripts to generate tcrBLOSUMs, PhysChemSim, and TopoSim matrices
  - *src/Dataset_assembly*: scripts to process TCR data and define blocks of epitope-specific TCRs
  - *src/Matrix_evaluation*: scripts to evaluate performance of tcrBLOSUMb when applied to clustering and alignment of distant TCRs

- *results/*: results and figures
  - */results/tcrBLOSUMmtx*: **tcrBLOSUM** matrices
    - *tcrBLOSUM_all_alpha.tsv* - tcrBLOSUMa, tcrBLOSUM matrix for CDR3 alpha
    - *tcrBLOSUM_all_beta.tsv* - tcrBLOSUMb, tcrBLOSUM matrix for CDR3 beta
  - */results/otherMTXs*: alternative amino acid substitution matrices
    - **PhysChemSimMtx.tsv** - a matrix which reflects the similarity between amino acids based on their physicochemical properties
    - **TopoSimMtx.tsv** - a matrix which reflects the similarity in atomic composition and atomic neighbourhoods between amino acids
