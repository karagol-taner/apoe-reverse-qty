# The Apolipoprotein E N-Terminal Bundle Rendered Membrane-Compatible by Reverse-QTY Conversion Without Loss of Fold

A substitution rule that reads nothing but secondary structure can rewrite 22 surface positions of a four-helix bundle, raise its apolar surface by 59%, and leave the fold intact.

This repository contains the source code and datasets for the computational design and analysis presented in the manuscript. It covers the reverse-QTY conversion, the AlphaFold3 structure comparison, the all-atom molecular dynamics analysis, and every figure in the paper.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Citation
If you use this framework in your research, please cite:

* Karagöl, T., & Karagöl, A. (2026). The Apolipoprotein E N-Terminal Bundle Rendered Membrane-Compatible by Reverse-QTY Conversion Without Loss of Fold.

## Usage

* **`codes/1-sequence_structure/`**: Sequence composition, hydropathy, hydrophobic moment, solvent accessibility and DSSP assignment from the AlphaFold3 models, and the independent re-derivation of that assignment from the two experimental structures.
* **`codes/2-molecular_dynamics/`**: Trajectory analysis for the three simulated systems: RMSD, radius of gyration, helix content, RMSF, solvent-accessible surface, lipid contacts, bilayer thickness and first-shell composition.
* **`codes/3-figures/`**: Figure generation. The scripts read only the JSON in `data/analysis_json/`, so the figures rebuild without re-running the analysis. `mkstruct.py` ray-traces the ribbons in Figure 2 and needs PyMOL.
* ! Note: the molecular dynamics analysis reads trajectories that are not in this repository. See **Data** below.
* **`data/AF3/`**: Both AlphaFold3 jobs, with all five models per sequence and their confidence files (pTM, pLDDT, PAE).
* **`data/Structures/`**: The experimental structures used for comparison, 1LPE (X-ray, 2.25 Å) and 2L7B (NMR, 20 conformers).
* **`data/analysis_json/`**: The analysis outputs every reported number is drawn from.
* **`data/MD_inputs/`**: Everything needed to rerun the three simulations: the CHARMM-GUI starting coordinates and topology for each system, and every GROMACS parameter file from minimisation through production.
* **`figures/`**: The seven published figures as vector PDF.
* **`rQTY_substitutions.csv`**: All 28 candidate glutamine, threonine and tyrosine positions in the segment, with DSSP code, helix, whether converted, and relative solvent accessibility in both sequences.

## Data

Three systems were simulated: the native segment in water (100 ns), the variant in a six-component neuronal bilayer (50 ns), and the variant in a highly mobile membrane mimetic (90 ns). The starting coordinates, topologies and every GROMACS parameter file are in `data/MD_inputs/`, so the simulations can be repeated from scratch. The quantities computed from the resulting trajectories are in `data/analysis_json/`, so every reported number and every figure can be checked without repeating them.

The production trajectories themselves are 337 MB and are not in this repository. They are available from the corresponding authors on request, and will be deposited in a public archive on publication.

## A note on residue numbering

The AlphaFold3 segment models are numbered 1 to 157 in their own coordinates, while 1LPE is numbered in mature ApoE coordinates. Every script derives the offset for each file by exact sequence match and asserts it over the whole overlap before measuring anything, because getting it wrong shifts every selection by ten residues with no visible error. Superposition RMSDs reported by PyMOL are re-derived independently by a Kabsch fit over explicitly matched residue numbers, and the scripts stop if the two disagree.

## Abstract

Apolipoprotein E (ApoE) performs most of its pathologically relevant biology at lipid interfaces, yet design efforts target its receptor binding, interdomain geometry or abundance rather than the lipid affinity of its N-terminal bundle. We asked whether the reverse-QTY (rQTY) code can make that bundle membrane-compatible without disturbing its fold. Within mature residues 11–167 we converted every glutamine, threonine and tyrosine lying inside an α-helix (Q→L, T→V, Y→F) and left those in turns and loops unchanged: 22 substitutions over 14.0% of the segment, sparing the receptor recognition region, moving the grand average of hydropathicity from −0.757 to +0.123. Across five AlphaFold3 models per sequence the variant bundle superimposes on the native at 0.81 Å mean Cα root-mean-square deviation with indistinguishable helix content, while apolar solvent-accessible surface rises 59% at constant total surface and burial. Hydrophobic moment increases in five of six helical segments, so amphipathicity is preserved. In all-atom molecular dynamics the native segment holds its fold in water over 100 ns; the variant stays folded and essentially fully helical in a six-component neuronal bilayer (50 ns) and a highly mobile membrane mimetic (90 ns), thinning the bilayer locally by 10.9 Å. The first lipid shell is only mildly biased in composition, and the two membrane models disagree. Each system was run once with the protein embedded at the build stage, so quantifying affinity will require matched native trajectories and partitioning free-energy calculations. rQTY thus acts as a geometry-preserving, residue-resolved control on the surface chemistry of a soluble helical bundle.

Keywords: Reverse-QTY code, Apolipoprotein E, Protein design, Amphipathic helix, Molecular dynamics, Alzheimer's disease
