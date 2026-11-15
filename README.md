# Salivary and plasma cytokine kinetics during experimental human endotoxaemia

Analysis code for the article *"Salivary and plasma cytokine kinetics during experimental human
endotoxaemia: macrophage inflammatory protein-1α and -1β as candidate non-invasive biomarkers"*
(Relouw FJA, Waalders NJB, Koolen J, Gerretsen J, Taal HR, Pickkers P, van Riel NAW, Kox M).
It makes the data figures and all tables of the article and its appendix from the study data.

## How to run

1. Put the two data files in `data/` (see `data/README.md`).
2. Install the Python packages once: `pip install -r python/requirements.txt`.
3. In MATLAB, open `run_all.m` and press **Run** (15–30 min).

By default the analysis starts from the model fit of the article (`fit/endotox_fit.mat`); all
other steps are recomputed. Set `REFIT = true` in `run_all.m` to refit the model with its 200
bootstrap refits (about 3 h with 16 parallel workers).

## Requirements

* MATLAB R2024b (recent releases should work) with the Optimization Toolbox. The Parallel
  Computing Toolbox is optional and speeds up the fits.
* Python 3.10 or later with the packages in `python/requirements.txt`.

## Files

| File | Content |
|---|---|
| `run_all.m` | runs everything |
| `data/README.md` | the two data files the code expects (not included) |
| `endotoxemia_model.m` | data import and cleaning, the compartmental model (ordinary differential equations), the three-stage fit and the bootstrap |
| `figures_tables.m` | figures, the shift analysis and the data and model tables |
| `hs.m` | figure style and export |
| `python/tables.py` | Table 1, values below the LLOQ, physicochemical properties (`sequences.fasta`: UniProtKB release 2026_04) |
| `fit/endotox_fit.mat` | model fit of the article: parameters, simulated curves and bootstrap results (no individual participant data) |

## Output

`results/figures/`: Figures 2–6 and S2–S8 (PDF, SVG, EPS, and 600-dpi PNG and TIFF).

`results/tables/`:

| File | Article |
|---|---|
| `table1_demographics.csv` | Table 1 |
| `shift_analysis.csv` | Table 2 |
| `kinetics.csv` | Table 3, Table S10 |
| `saliva_below_lloq.csv`, `below_lloq_by_time.csv` | Table S1 |
| `below_lloq_overall.csv` | Methods, values below the LLOQ |
| `model_metrics.csv` | Table S2 (all derived model quantities with 95% CIs) |
| `parameters.csv`, `baselines.csv` | Tables S3–S4 |
| `spread.csv` | Table S6 |
| `robustness.csv` | Tables S7 and S11 |
| `lags_leave_one_group_out.csv` | Table S7 |
| `peaks.csv` | Tables S8 and S10 |
| `half_lives.csv` | Table S9 |
| `physicochemical.csv` | Tables S12–S13 |
| `arms_summary.csv` | Results, per treatment group |

## Notes

* Random numbers use fixed seeds (2026), so the results are reproducible. A refit on another
  machine or MATLAB release can differ in the last digits of the optimiser output.
* Only the study day is analysed: the LPS challenge and its 24-h sample (coded D0 and D1 in the
  data file).

## Data and licence

The participant data are not included; see the data sharing statement of the article.
Code licence: MIT (see `LICENSE`).
