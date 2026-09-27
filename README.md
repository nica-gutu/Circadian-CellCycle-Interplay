# Circadian clock–cell cycle interplay

Thesis-specific analysis and modeling repository accompanying the PhD thesis:

**Nica Gutu — _Deciphering chronobiological regulation of cell proliferation and drug responses: Insights from the circadian clock and p53-p21 dynamics_**

Thesis: https://edoc.hu-berlin.de/items/8808e9d0-e627-4f01-83fc-6f2507daf420

## Purpose of this repository

This repository contains a **thesis-specific reanalysis and reconstruction** of the circadian clock–cell cycle project.

The underlying biological project overlaps with the work published as:

**Gutu, N. et al. “Circadian coupling orchestrates cell growth.” _Nature Physics_ 21, 768–777 (2025).**  
https://doi.org/10.1038/s41567-025-02838-4

However, this repository is **not simply a copy of the publication analysis repository**. Analyses and figures were independently recreated for the thesis so that the thesis could present the work without directly reusing the published figure assets. It also contains additional analyses and visualizations that were useful for the thesis but do not appear in the final paper.

For the code corresponding directly to the published manuscript figures, see:

https://github.com/nica-gutu/Circadian-clock-cell-cycle

## Scientific focus

The project investigates how **extracellular circadian coupling** influences synchronization among individual cellular clocks, coordination between the circadian clock and the cell cycle, and proliferation at single-cell and population scales.

The thesis analysis emphasizes the quantitative relationship between:

- circadian synchronization across cells;
- intracellular circadian–cell-cycle phase coordination;
- single-cell division timing;
- population-level growth;
- mathematical descriptions of coupled oscillatory systems.

## Repository structure

### `Circadian-CellCycle/`

Thesis-specific circadian and clock–cell-cycle analyses, including:

- circadian amplitude and period;
- fraction of oscillatory cells;
- phase coherence and phase locking;
- extracellular versus intracellular coupling;
- Poincaré-oscillator simulations;
- entrainment regions and parameter sweeps;
- experimental/model comparisons;
- additional phase-coherence and spatially resolved analyses.

### `Population-proliferation/`

Population-level growth analyses, including:

- growth under altered circadian coordination;
- density-dependent effects;
- comparison of wild-type and clock-disrupted conditions.

### `Single-cell-proliferation/`

Single-cell proliferation analyses, including:

- intermitotic-time distributions;
- oscillatory-period distributions;
- proliferation trajectories and related visualizations.

## Computational approaches

The repository contains analyses based on:

- coupled nonlinear oscillator modeling;
- numerical simulation of circadian and cell-cycle dynamics;
- circadian period, amplitude, phase, and phase-coherence quantification;
- entrainment and parameter-space analysis;
- time-frequency analysis using `pyBOAT`;
- single-cell and population-level proliferation measurements;
- comparison of experimental observations with mechanistic models.

## Relationship to the publication repository

The two repositories serve different purposes:

- **`Circadian-clock-cell-cycle`** — code corresponding to the published _Nature Physics_ study and its figures.
- **`Circadian-CellCycle-Interplay`** — thesis-oriented reconstruction, alternative visualizations, and additional analyses used to present and extend the project in the PhD thesis.

Because some analyses address the same scientific questions, overlap in computational methods and script logic is expected.

## Software

The analyses are implemented in Python and use packages including:

```text
numpy
pandas
matplotlib
scipy
pyboat
```

Some scripts may use additional scientific-Python packages depending on the analysis.

## Data and reproducibility

This repository primarily contains thesis analysis and figure-generation code. The associated experimental data for the published circadian-coupling study are available on Figshare:

https://doi.org/10.6084/m9.figshare.28375358.v1

Some scripts retain paths from the original research environment and require the corresponding source data and local path configuration before execution.

## Related publication

> Gutu, N., Nordentoft, M.S., Kuhn, M. et al. **Circadian coupling orchestrates cell growth.** _Nature Physics_ 21, 768–777 (2025). https://doi.org/10.1038/s41567-025-02838-4

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
