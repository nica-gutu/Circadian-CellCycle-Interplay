# Circadian coupling, cell-cycle coordination, and cell growth

Analysis and modeling code related to:

**Gutu, N. et al. “Circadian coupling orchestrates cell growth.” _Nature Physics_ 21, 768–777 (2025).**  
https://doi.org/10.1038/s41567-025-02838-4

## Overview

This project investigates how **extracellular coupling between single-cell circadian oscillators** influences coordination between the circadian clock, the cell cycle, and collective cell growth.

The study combines coupled-oscillator theory with long-term population and single-cell recordings, live-cell imaging, perturbation experiments, and machine-learning-assisted image analysis. The central question is whether tissue-level circadian synchronization merely synchronizes clocks or also changes intracellular clock–cell-cycle coordination and, consequently, proliferation.

## Experimental and computational framework

The study uses human U2OS cell systems with circadian and cell-cycle reporters. Single-cell recordings quantify circadian dynamics and division timing, while population-level measurements quantify confluence and growth.

Circadian coupling is perturbed in complementary ways, including:

- pharmacological inhibition of TGF-β signaling with **LY2109761**;
- changes in cell density;
- genetic disruption of the circadian clock, including **Cry1/Cry2 double-knockout** cells.

A coupled-oscillator model is used to separate **extracellular circadian coupling** from **intracellular coupling between the circadian clock and cell cycle**, and to study how both determine synchronization and phase locking.

Single-cell image sequences were processed using an automated pipeline that included supervised pixel/object classification, segmentation, tracking, and signal extraction with **ilastik**.

## Main findings represented by the analyses

- Increasing extracellular coupling increases synchronization among cellular circadian oscillators.
- Weakening extracellular circadian coupling accelerates circadian desynchronization.
- Loss of circadian synchronization disrupts phase coordination between the circadian clock and cell cycle within individual cells.
- Reduced circadian coupling is associated with altered cell-cycle timing and impaired collective growth.
- Coherent circadian populations display **oscillatory growth dynamics**, linking circadian synchronization to tissue-level proliferation.
- Genetic disruption of the core clock removes the coupling-dependent growth effect, supporting a specific role for circadian-clock regulation rather than a nonspecific drug effect.

Together, the results support a multiscale mechanism in which cell-to-cell circadian communication influences intracellular clock–cell-cycle coordination and thereby regulates population growth.

## Repository structure

### `Circadian-CellCycle/`

Circadian-clock and cell-cycle analyses, including:

- circadian period and amplitude;
- phase coherence and phase locking;
- coupled-oscillator simulations;
- extracellular/intracellular coupling parameter sweeps;
- entrainment regions;
- comparison of model and experimental dynamics.

### `Population-proliferation/`

Population-level proliferation and growth analyses, including the effects of altered circadian coordination and comparison of wild-type and clock-disrupted conditions.

### `Single-cell-proliferation/`

Single-cell proliferation analyses, including intermitotic-time distributions, oscillatory periods, and division-related measurements.

## Computational methods

Methods represented in the repository include:

- coupled nonlinear oscillator models;
- numerical simulation of interacting circadian and cell-cycle oscillators;
- parameter-space and entrainment-region analysis;
- circadian amplitude, period, phase, and phase-coherence quantification;
- time-frequency analysis using `pyBOAT`;
- single-cell and population-level proliferation analysis;
- fitting of population growth curves;
- integration of experimental recordings with mathematical modeling.

## Software

The published study reports Python 3.8.8 and the following core analysis packages:

```text
numpy 1.23.3
pandas 1.4.4
matplotlib 3.6.3
seaborn 0.11.2
scipy 1.9.1
scikit-image 0.20.0
pyBOAT 0.9.1
```

The imaging workflow additionally used **ilastik** for supervised segmentation, tracking, and signal extraction.

## Data and reproducibility

The experimental raw data and processed tables associated with the study are publicly available on Figshare:

**Dataset:** https://doi.org/10.6084/m9.figshare.28375358.v1

The publication lists the primary study code at:

- https://github.com/Granada-Lab/Circadian-clock-cell-cycle
- https://github.com/Granada-Lab/Automatic-single-cell-tracking-of-fading-objects

Some scripts in this personal repository retain analysis-specific paths from the original research environment and may require path adjustment before execution.

## Citation

> Gutu, N., Nordentoft, M.S., Kuhn, M. et al. **Circadian coupling orchestrates cell growth.** _Nature Physics_ 21, 768–777 (2025). https://doi.org/10.1038/s41567-025-02838-4

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
