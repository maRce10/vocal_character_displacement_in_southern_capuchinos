# Vocal character displacement in southern capuchinos

[![Website](https://img.shields.io/badge/website-analysis%20reports-408AB4)](https://marce10.github.io/vocal-character-displacement-capuchinos/)
[![Publish website](https://github.com/maRce10/vocal-character-displacement-capuchinos/actions/workflows/publish.yml/badge.svg)](https://github.com/maRce10/vocal-character-displacement-capuchinos/actions/workflows/publish.yml)

Code and data for the analysis of vocal character displacement in southern capuchino seedeaters (*Sporophila* spp.). We evaluate whether songs of different species are more or less divergent in sympatry than in allopatry, how much each species' song varies across its range, and which species' populations shifted where they live with congeners. Simple and complex songs are analyzed separately.

## Website

The rendered analysis reports, with all code, results and figures, are available at **[marce10.github.io/vocal-character-displacement-capuchinos](https://marce10.github.io/vocal-character-displacement-capuchinos/)**:

- [Data preparation](https://marce10.github.io/vocal-character-displacement-capuchinos/scripts/0_data_preparation.html): PCA of element- and song-level acoustic features and pairwise acoustic and geographic distances
- [Step 1: Population structure](https://marce10.github.io/vocal-character-displacement-capuchinos/scripts/1_population_structure.html): within- vs. between-population song variation and isolation by distance
- [Step 2: Sympatry](https://marce10.github.io/vocal-character-displacement-capuchinos/scripts/2_sympatry_and_heterospecific_distances.html): acoustic distance between heterospecific songs in sympatry vs. allopatry
- [Step 3: Population shifts](https://marce10.github.io/vocal-character-displacement-capuchinos/scripts/3_population_shifts.html): whether populations that co-occur with a congener are shifted towards or away from its song

## Repository structure

```
_quarto.yml, index.qmd, styles.css            Quarto website configuration (cached code results in _freeze/)
.github/workflows/publish.yml                 GitHub Actions workflow that builds the website and deploys it to GitHub Pages
scripts/
  0_data_preparation.qmd                      PCA and pairwise distances; saves data/processed/step_inputs_*.rds
  1_population_structure.qmd                  Step 1 models
  2_sympatry_and_heterospecific_distances.qmd Step 2 models (sympatry-only, PCA-based and per-feature)
  3_population_shifts.qmd                     Step 3 models
  functions.R                                 Custom functions sourced by every report
  qmd.css                                     Style sheet for the reports
data/
  raw/                                        Song- and element-level acoustic measurements and individual coordinates
  processed/
    step_inputs_simple.rds, step_inputs_complex.rds
                                              Pairwise data sets used by Steps 1-3 (created by 0_data_preparation.qmd)
    fits/                                     Fitted brms models (not on GitHub, see below)
```

## Rendering the website

Code results are cached in `_freeze/` (`execute: freeze: auto`), so the GitHub Actions workflow only assembles the website and does not need R. After editing a report, render it locally, from the project root:

```
quarto render scripts/1_population_structure.qmd
```

and commit the updated `_freeze/` folder together with the report. Every push to `master` rebuilds and redeploys the website.

Run `0_data_preparation.qmd` before the step reports, as it creates the data they read.

## Files not hosted on GitHub

- `data/processed/fits/`: fitted brms models (several hundred MB). They are recreated when the reports are rendered.
- `manuscript/` and `archive/` are kept locally only.
