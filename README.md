# Integrative Ethnopharmacological and In-Silico Strategy for Identifying Indonesian Jamu to Enhance Athletic Stamina

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.13](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/downloads/)
[![AutoDock Vina](https://img.shields.io/badge/AutoDock%20Vina-1.2.5-green.svg)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![ChimeraX](https://img.shields.io/badge/UCSF%20ChimeraX-1.11-orange.svg)](https://www.rbvi.ucsf.edu/chimerax/)

## Overview

This repository contains the complete computational pipeline, data, and results for the molecular docking study:

> **"Integrative Ethnopharmacological and In-Silico Strategy for Identifying Indonesian Jamu to Enhance Athletic Stamina"**

### Summary

We performed systematic molecular docking of **31 bioactive phytochemicals** against **4 human stamina-associated proteins** (ANDR_HUMAN, ESR1_HUMAN, HIF1A_HUMAN, NFKB1_HUMAN) using AutoDock Vina 1.2.5. Protein structures were obtained from AlphaFold2, prepared using UCSF ChimeraX, and validated through native ligand re-docking. Binding interactions were analyzed using ChimeraX and LigPlot+.

---

## Repository Structure

```
molecular-docking-phytochemicals/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── requirements.txt             # Python dependencies
├── environment.yml              # Conda environment
│
├── scripts/
│   ├── 01_protein_prep.cxc      # ChimeraX protein preparation script
│   ├── 02_run_docking.py        # Main docking pipeline (all 31 compounds)
│   ├── 03_parse_results.py      # Parse and summarize Vina outputs
│   ├── 04_generate_complex.py   # Generate receptor-ligand complex PDB for LigPlot+
│   ├── 05_publication_viz.py    # Generate ChimeraX .cxc visualization scripts
│   └── utils.py                 # Shared utility functions
│
├── config/
│   └── batch_config.csv         # All 31 compounds with SMILES, targets, grid params
│
├── data/
│   ├── compounds_final.csv      # Full compound dataset with ADMET properties
│   └── pocket_coordinates.json  # Validated binding site coordinates per protein
│
├── results/
│   ├── best_poses_summary.csv   # Best docking pose per compound-protein pair
│   ├── combined_docking_results.csv  # All 9 poses per pair (279 rows)
│   └── interactions/
│       └── hbond_contacts_all.csv    # H-bond and contact data for all pairs
│
└── docs/
    ├── methods.md               # Detailed methods (plain text)
    └── CITATION.cff             # Citation file
```

---

## Target Proteins

| Protein | Full Name | UniProt | PDB (Validation) | Native Ligand |
|---------|-----------|---------|------------------|---------------|
| ANDR_HUMAN | Androgen Receptor | P10275 | — | Testosterone (CID: 6013) |
| ESR1_HUMAN | Estrogen Receptor α | P03372 | 1GWR | Estradiol (CID: 5757) |
| HIF1A_HUMAN | HIF-1α | Q16665 | — | 1,4-DPCA (CID: 459803) |
| NFKB1_HUMAN | NF-kB p50 | P19838 | 8TQD | JMR inhibitor |

---

## Key Results

| Rank | Compound | Target | Affinity (kcal/mol) | vs Native |
|------|----------|--------|---------------------|-----------|
| 1 | Oestrone | ESR1 | −9.45 | +1.29 |
| 2 | Equol | ESR1 | −8.35 | +0.19 |
| 3 | Estriol | ESR1 | −8.16 | = native |
| 4 | Chrysin | ESR1 | −8.03 | — |
| 5 | Withanolide | NFKB1 | −7.66 | +3.24 |
| 6 | Estradiol | HIF1A | −7.39 | +1.18 |

---

## Installation & Usage

### Prerequisites

- Python ≥ 3.10
- [UCSF ChimeraX 1.11](https://www.rbvi.ucsf.edu/chimerax/download.html)
- [AutoDock Vina 1.2.5](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
- Conda (recommended) or pip

### 1. Clone repository

```bash
git clone https://github.com/[your-username]/molecular-docking-phytochemicals.git
cd molecular-docking-phytochemicals
```

### 2. Create conda environment

```bash
conda env create -f environment.yml
conda activate docking
```

Or with pip:

```bash
pip install -r requirements.txt
```

### 3. Install AutoDock Vina

**macOS (Apple Silicon):**
```bash
# Download arm64 binary
curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_arm64 -o vina
chmod +x vina
sudo mv vina /usr/local/bin/vina
```

**macOS (Intel):**
```bash
curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_x86_64 -o vina
chmod +x vina && sudo mv vina /usr/local/bin/vina
```

**Linux:**
```bash
sudo apt install autodock-vina
```

### 4. Prepare protein structures

Download AlphaFold structures for the 4 target proteins:

```bash
# ANDR_HUMAN (UniProt: P10275)
wget https://alphafold.ebi.ac.uk/files/AF-P10275-F1-model_v4.pdb -O protein/ANDR_HUMAN.pdb

# ESR1_HUMAN (UniProt: P03372)
wget https://alphafold.ebi.ac.uk/files/AF-P03372-F1-model_v4.pdb -O protein/ESR1_HUMAN.pdb

# HIF1A_HUMAN (UniProt: Q16665)
wget https://alphafold.ebi.ac.uk/files/AF-Q16665-F1-model_v4.pdb -O protein/HIF1A_HUMAN.pdb

# NFKB1_HUMAN (UniProt: P19838)
wget https://alphafold.ebi.ac.uk/files/AF-P19838-F1-model_v4.pdb -O protein/NFKB1_HUMAN.pdb
```

### 5. Prepare receptors in ChimeraX

Open ChimeraX and run for each protein:
```
open protein/ANDR_HUMAN.pdb
delete @@bfactor<70
delete solvent
delete ligand
addh hbond true
addcharge method am1-bcc
save results/ANDR_HUMAN_prepped.pdb
close #1
```
*(Repeat for ESR1, HIF1A, NFKB1)*

### 6. Run docking pipeline

```bash
python scripts/02_run_docking.py --config config/batch_config.csv --out results/
```

### 7. Parse and summarize results

```bash
python scripts/03_parse_results.py --results results/ --out results/best_poses_summary.csv
```

---

## Data Availability

- **Protein structures**: AlphaFold2 (https://alphafold.ebi.ac.uk)
- **Compound SMILES**: PubChem (https://pubchem.ncbi.nlm.nih.gov)
- **Bioactivity data**: ChEMBL (https://www.ebi.ac.uk/chembl)
- **Phytochemical data**: KNApSAcK (http://www.knapsackfamily.com)
- **Docking results**: Available in `results/` directory of this repository

> **Note**: Due to file size constraints, raw PDBQT docking pose files are not included. They can be reproduced by following the pipeline above. The `results/best_poses_summary.csv` contains all binding affinity scores.

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{[AuthorLastName][Year],
  title   = {Molecular Docking Study of Bioactive Compounds from Traditional Medicinal Plants Against Cancer-Associated Target Proteins},
  author  = {[Author Names]},
  journal = {[Journal Name]},
  year    = {[Year]},
  doi     = {[DOI]}
}
```

---

## Software Citations

Please also cite the tools used:

- **AutoDock Vina**: Eberhardt, J. et al. (2021). *J. Chem. Inf. Model.*, 61(8), 3891–3898. https://doi.org/10.1021/acs.jcim.1c00203
- **UCSF ChimeraX**: Meng, E.C. et al. (2023). *Protein Sci.*, 32(11), e4792. https://doi.org/10.1002/pro.4792
- **AlphaFold2**: Jumper, J. et al. (2021). *Nature*, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2
- **RDKit**: Landrum, G. (2024). RDKit: Open-Source Cheminformatics. https://www.rdkit.org
- **meeko**: Forli Lab (2023). meeko: Preparation of small molecules for AutoDock. https://github.com/forlilab/meeko
- **LigPlot+**: Laskowski, R.A. & Swindells, M.B. (2011). *J. Chem. Inf. Model.*, 51(10), 2778–2786.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

Matthew Valentino Tambunan — tambunan.matthewv@gmail.com
Department of Pharmacy, Faculty of Mathematics and Natural Sciences, Udayana University, Bali, Indonesia
