# A Space–Time Dual-Domain Test of the Source-Superposition Law in Volcanic Systems

**Repository:** https://github.com/Ruqing1963/volcano-spacetime-rmt

This repository accompanies the paper:

**"A Space–Time Dual-Domain Test of the Source-Superposition Law in Volcanic
Systems: Temporal Recharge Repulsion, Spatial Supply-Shadow Repulsion, and the
Spectral Superposition Theorem Across Four Natural Targets"**
by Ruqing Chen (GUT Geoservice Inc., Montreal).

## Summary

This is the capstone of a research program applying Random Matrix Theory (RMT)
to Earth-system rhythms. The program's unifying principle:

> **A single isolated long-memory process yields GOE/GUE level repulsion; the
> superposition of many independent sources yields Poisson/clustered statistics.**

Earlier papers established this law separately in the **temporal** domain
(tectonic rhythms, seismicity, mantle-plume and metallogenic event timing) and
in the **spatial** domain (cyclostratigraphic bed spacing; orogenic-gold vein
clustering). This paper closes the program by probing a single volcanic *class*
in **both** domains at once, across four independent natural targets, using only
the parameter-free spacing-ratio statistic ⟨r⟩ and the coefficient of variation
CV — with no tunable preprocessing.

The four targets realize all four logical cases of the law:

- **(A) Time, single source — Mount Etna.** 86 flank eruptions (1169–2021 CE).
  An isolated magma chamber recharges between eruptions → temporal repulsion.
- **(B) Time, many sources — six-volcano pool.** Per-source-unfolded intervals
  from six independent volcanoes, randomly interleaved → spectral superposition
  theorem → relaxation to Poisson.
- **(C) Space, single arc — Cascades.** 17 stratovolcanoes (Baker→Lassen).
  Magma-supply shadows exclude neighbours → spatial GOE repulsion.
- **(D) Space, complex multi-rift — Iceland.** 67 vents from seven rift systems,
  projected unfiltered onto one transect. Structural inheritance → clustering.

## Key results

| Target | Domain | Type | Dataset | N | ⟨r⟩ | β̂ | CV | Best fit | Pass |
|---|---|---|---|---|---|---|---|---|---|
| A | Time | Single | Etna flank | 85 | 0.518 | 0.88 | 1.16 | ~GOE | ✓ |
| B | Time | Multi | 6-volcano pool | 455 | 0.470 | 0.56 | 1.27 | Poisson | ✓ |
| C | Space | Arc | Cascade arc | 16 | 0.594 | 1.86 | 0.56 | GOE | ✓ |
| D | Space | Complex | Iceland rifts | 66 | 0.672 | 2.69 | 1.59 | Clustered | ✓ |

Reference values: Poisson ⟨r⟩ = 0.386, CV = 1.0; GOE ⟨r⟩ = 0.536, CV ≈ 0.52;
GUE ⟨r⟩ = 0.603.

Highlights:

- **Target B** drops ⟨r⟩ from a per-source mean of **0.574** to **0.470** once
  the six sources are interleaved — a direct confirmation of the superposition
  theorem.
- **Target C** is the cleanest single result: CV = 0.56 sits almost exactly on
  the GOE reference, and the GOE Kolmogorov–Smirnov fit **cannot be rejected**
  (p = 0.79), with ~54 km quasi-regular along-arc spacing.
- **Targets A + D** reproduce, in volcanic form, the **time–space asymmetry**
  first reported for orogenic gold: temporally repulsive, spatially clustered.

## Methodological notes

- **Per-source unfold then pool (Target B).** Each volcano's intervals are
  normalized by its own mean rate, then pooled and **randomly shuffled** to
  destroy within-source ordering. Pooling without shuffling, or naive
  chronological merging of raw dates, fabricates artefactually high ⟨r⟩ (0.62
  and 0.78 respectively); only the decorrelated result is reported.
- **No declustering (Targets C, D).** Minimum-separation declustering fabricates
  spurious GOE statistics from random points and is never used. All spatial
  statistics are computed on raw projected positions.
- **Read ⟨r⟩, β̂ and CV together.** Target D shows why: elevated ⟨r⟩ = 0.672
  *and* CV ≫ 1 jointly diagnose dual-scale structure (regular within clusters,
  clustered between them); either statistic alone would mislead.

## Repository layout

```
code/      volcano_rmt_pipeline.py   reproducible four-target analysis + figure
           export_data.py            writes the data/ CSVs
data/      target_A_etna_flank_eruptions.csv      86 Etna eruption years
           target_B_multivolcano_eruptions.csv    461 eruptions, 6 volcanoes
           target_C_cascade_arc_positions.csv     17 Cascade volcano positions
           target_D_iceland_vents.csv             67 Iceland vents + systems
figures/   volcano_spacetime_rmt.pdf  four-panel result figure (vector)
           volcano_spacetime_rmt.png  raster version
paper/     volcano_spacetime_rmt.tex  LaTeX source
           volcano_spacetime_rmt.pdf  compiled paper
results/   volcano_rmt_results.json   full numerical results (all four targets)
```

## Data sources

Primary eruption chronologies derive from the **Smithsonian Global Volcanism
Program** (Holocene eruption database) supplemented by published
volcano-specific records (e.g. Branca & Del Carlo 2005 for Etna). Cascade-arc
positions are projected from GVP coordinates onto the local arc strike
(N10°W); Icelandic vents are projected onto a single N25°E mega-transect. The
processed datasets used in the analysis are provided in `data/`.

## Reproduce

```bash
pip install -r requirements.txt
cd code && python3 volcano_rmt_pipeline.py
```

This regenerates `figures/volcano_spacetime_rmt.{png,pdf}` and
`results/volcano_rmt_results.json`. To rebuild the CSV datasets, run
`python3 export_data.py` from `code/`.

## Related work in this program

- Temporal racetracks (faults, mantle plumes, metallogeny): GOE level repulsion
  from single long-memory sources.
- Spatial Part I (cyclostratigraphy, GUE/GOE repulsion) and Part II (Sigma
  orogenic-gold veins, spatial clustering; the time–space asymmetry).

## License

MIT (see LICENSE).
