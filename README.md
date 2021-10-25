# Mediators_of_HeadMotion
This repository contains the scripts for the analyses of mediators between body mass index (BMI) and Head Motion (HM).

## Main Analysis

`extract_physio_correlates.m` is a preparatory script to extract HM, Bold signal fluctuations (DVARS) and physiological parameters (respiratory volume per time, respiratory rate, heart rate) from the files created during preprocessing.
See [here](https://github.com/fBeyer89/life_followup_preproc) for the preprocessing code which generated `rest_realigned.nii.gz.par` (for HM), `rest2anat_dvars` (for DVARS) and `*.mat` (for physiological parameters).

The main analyses steps (loading the data, preparing the data and running mediation models ) is in `Mediation_BMI_HeadMotion.Rmd`.
The Rmd file was rendered to `Mediation_BMI_HeadMotion.html`.

The `Rmd` uses `create_mediation_model.R`  and  `create_regression_model.R`  to build the respective models and run them for the different mediators. In the `Rmd`, all mediation plots and the flowchart of subject inclusion/exclusion are created (`flowchart.html`, `graf_revised.png`, `Mediators_of_HeadMotion/Flowchart.png`).

### Additional scripts not actively used in manuscript

`plot_correlations_for_single_subject.m` allows to plot correlations of HM, DVARS and physiological parameters from the PhysIO toolbox per subject.

`compare_mediation_model.R` was a script to compare another formulation of mediation with the one finally used in manuscript.

`BIS_preparation` contains an old `.sps` script which informed our analysis of the BIS.
