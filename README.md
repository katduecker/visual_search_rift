# visual_search_rift

Code for the experimental paradigm, MEG/EEG analysis pipelines, and statistics behind two studies that combine Magnetoencephalography (MEG) and Rapid Invisible Frequency Tagging (RIFT) in a classic visual search task. Both studies analyse the same dataset of 31 participants searching for a single target "T" among 16 or 32 "L" distractors, under guided and unguided conditions.

## Associated papers

**Paper 1 (published).** Duecker, K., Shapiro, K. L., Hanslmayr, S., Griffiths, B. J., Pan, Y., Wolfe, J. M., & Jensen, O. (2025). Guided visual search is associated with target boosting and distractor suppression in early visual cortex. *Communications Biology*, 8, 912. https://doi.org/10.1038/s42003-025-08321-3

This study shows that knowing the target colour modulates neuronal excitability in early visual cortex in line with a priority map. RIFT responses to the target colour are boosted and responses to the distractor colour are suppressed in guided relative to unguided search. The effect is established with two complementary measures: magnitude-squared coherence (`coherence/`) and a single-trial General Linear Model of the RIFT response (`glm_spec/`).

**Paper 2 (in revision).** Duecker, K., Shapiro, K. L., Hanslmayr, S., Griffiths, B. J., Quinn, A. J., Wolfe, J. M., Pan, Y., Pastuszak, A., & Jensen, O. Higher baseline alpha power is associated with faster responses in visual search.

This study re-analyses the same data to ask whether occipital alpha oscillations support search. Higher alpha power before and during the search correlates with faster responses, and this holds after controlling for time-on-task using a GLM-spectrum approach (`glm_spec/`, `alpha/`). The increase in alpha power does not correlate with a reduced RIFT response, so the analyses provide no direct evidence for feature-specific gain modulation by alpha.

## Data availability

Raw MEG data are archived on Zenodo: https://zenodo.org/records/17752242

## Repository layout

The folders follow the order of the processing pipeline. Within each folder, script prefixes (`a_`, `b1_`, `c1_`, and so on) indicate the run order. Most analysis scripts take a subject index as input and are paired with a `.sh` wrapper for batch submission on an HPC cluster.

| Folder | Language | Contents |
|--------|----------|----------|
| `experiment/` | MATLAB (Psychtoolbox 3) | The visual search paradigm. Entry point is `a_exp.m`. Includes RIFT stimulation, EyeLink setup, and MEG triggering. |
| `mne_scripts/` | Python (MNE) | Noisy-sensor detection and Signal Space Separation ("Maxfilter"). Run before the MATLAB pipeline. |
| `matlab_scripts/preproc_meg/` | MATLAB (FieldTrip) | Trial definition, EDF/FIF/response-file merging, semi-automatic artefact rejection, ICA, photodiode-to-sine replacement, and sensor-of-interest selection. |
| `matlab_scripts/coherence/` | MATLAB (FieldTrip) | Magnitude-squared coherence pipeline for the RIFT response, per condition and at the single-trial level. (Paper 1 and 2) |
| `matlab_scripts/glm_spec/` | MATLAB | GLM analyses. Single-trial GLM of the RIFT response (Paper 1), the GLM-spectrum linking alpha power to reaction time with interaction terms (Paper 2), and the GLM linking the RIFT response to alpha power (Paper 2). |
| `matlab_scripts/alpha/` | MATLAB (FieldTrip) | Alpha-band time-frequency analysis, individual alpha frequency and sensor-of-interest identification, condition contrasts, and fast-versus-slow control analyses (Paper 2). |
| `matlab_scripts/source/` | MATLAB (FieldTrip) | T1-to-head alignment, forward model, and DICS beamformer source localisation of the RIFT and alpha responses. (Paper 1)|
| `matlab_scripts/gaze/` | MATLAB | Eye-movement and gaze-bias control analyses from the EyeLink data. (Paper 1)|
| `matlab_scripts/behavior/` | MATLAB | Signal detection (d') and reaction time summaries by condition. (Paper 1)|
| `matlab_scripts/cbrewer/` | MATLAB | ColorBrewer colormap utility (third party). |
| `r_scripts/` | R | Behavioural linear mixed models, eye-movement statistics. (Paper 1)|

## Pipeline order

1. `experiment/` Code to run RIFT experiment.
2. `mne_scripts/` cleans sensors and applies SSS.
3. `matlab_scripts/preproc_meg/` defines trials, rejects artefacts, runs ICA, and selects sensors of interest.
4. The analysis folders then run independently, depending on the question:
   - RIFT priority-map effects: `coherence/` and `glm_spec/` (Paper 1).
   - Alpha and reaction time: `alpha/` and `glm_spec/` (Paper 2).
   - Source localisation: `source/` (Paper 1).
   - Controls: `gaze/`, `behavior/` (Paper 1).
5. `r_scripts/` produces the behavioural and group-level statistics.

## Software requirements

- MATLAB R2019b.
- [FieldTrip](https://www.fieldtriptoolbox.org/) for MEG preprocessing, frequency analysis, source modelling, and cluster-based permutation tests.
- [Psychtoolbox 3](http://psychtoolbox.org/) for stimulus presentation.
- [MNE-Python](https://mne.tools/) for noisy-sensor detection and Maxfilter/SSS.
- R (the behavioural analyses were run in RStudio 1.1.456 with R 3.6.1).
- The GLM-spectrum analyses follow Quinn et al. (2024), *Imaging Neuroscience*, https://doi.org/10.1162/imag_a_00082.

## Notes on running the code

The scripts contain absolute paths from the Birmingham HPC environment (for example `/rds/projects/j/jenseno-visual-search-rft/`). Update the `pth` variables and the FieldTrip path at the top of each script to match your own setup. The `.sh` files are batch wrappers and assume a SLURM scheduler; adapt or ignore them if you run the scripts interactively.

## Citation

If you use this code, please cite the relevant paper:

```
Duecker, K., Shapiro, K. L., Hanslmayr, S., Griffiths, B. J., Pan, Y., Wolfe, J. M., & Jensen, O. (2025).
Guided visual search is associated with target boosting and distractor suppression in early visual cortex.
Communications Biology, 8, 912. https://doi.org/10.1038/s42003-025-08321-3
```

## Contact

Katharina Duecker, Department of Neuroscience, Brown University. katharina.duecker@gmail.com
