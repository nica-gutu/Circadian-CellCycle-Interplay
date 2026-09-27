# Circadian coupling, cell-cycle coordination, and cell growth

Analysis and modeling code supporting the study **“Circadian coupling orchestrates cell growth”** (*Nature Physics*, 2025).

**Publication:** https://doi.org/10.1038/s41567-025-02838-4

## Scientific question

This project investigates how communication between single-cell circadian oscillators influences coordination between the **circadian clock**, the **cell cycle**, and collective cell growth.

The study combines coupled-oscillator modeling with long-term live-cell measurements, perturbation experiments, quantitative time-series analysis, and machine-learning-assisted bioimage analysis. The central goal is to understand how loss of extracellular circadian synchronization changes clock–cell-cycle coordination within individual cells and how these changes propagate to tissue-level growth dynamics.

## Repository contents

The code is grouped by the biological scale or analysis represented in the study.

### `Circadian-CellCycle/`

Analysis and modeling of circadian-clock and cell-cycle dynamics, including:

- circadian period and amplitude measurements;
- phase coherence and phase-locking analyses;
- coupled-oscillator simulations;
- entrainment regions and parameter-space exploration;
- comparison of experimental and model-derived circadian properties.

### `Population-proliferation/`

Population-level proliferation analyses, including growth under altered circadian coordination and comparison of wild-type and clock-perturbed conditions.

### `Single-cell-proliferation/`

Single-cell proliferation analyses, including intermitotic-time distributions, oscillatory periods, and proliferation-related measurements.

The repository contains research-analysis scripts corresponding to analyses and figures from the study rather than a packaged software library.

## Computational approaches

Methods represented in the repository include:

- coupled nonlinear oscillator models;
- numerical simulation of interacting circadian and cell-cycle oscillators;
- single-cell and population-level time-series analysis;
- circadian amplitude, period, and phase-coherence quantification;
- wavelet-based signal analysis using `pyBOAT`;
- proliferation and intermitotic-time analysis;
- integration of experimental measurements with mathematical models.

The experimental workflow associated with the study also used machine-learning-based bioimage analysis with **ilastik** for single-cell tracking and feature extraction.

## Requirements

The scripts are written in Python. Core dependencies used across the repository include:

```text
numpy
pandas
matplotlib
pyboat
```

Some scripts may use additional scientific-Python packages depending on the analysis.

## Data and reproducibility

The repository primarily contains the analysis and modeling code. Raw experimental imaging and measurement datasets are **not included here**. Some scripts contain analysis-specific input/output paths and therefore require the relevant source data and local path configuration before they can be rerun.

The folder structure follows the main analysis components of the paper, making it possible to locate the code associated with circadian–cell-cycle dynamics, population growth, and single-cell proliferation separately.

## Citation

If you use this code or build on these analyses, please cite:

> Gutu, N., Nordentoft, M.S., Kuhn, M. et al. **Circadian coupling orchestrates cell growth.** *Nature Physics* 21, 768–777 (2025). https://doi.org/10.1038/s41567-025-02838-4

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
