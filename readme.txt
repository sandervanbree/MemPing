These files were used to run and analyse the long-term memory ping study by Sander van Bree, Abbie Sarah Mackenzie, and Maria Wimber.
Paper: [link pending] (open-access CC BY)

- The Analyse folder contains sequential analysis scripts (with a separate set of revisions analyses)
- The Dependencies folder contains scripts and files you need to run the analyses.
- The Experiment folder contains Psychopy material used for running the study
- The Simulation folder contains material for running analysis-validating simulations
- A Data folder needs to be created locally; the EEG and behavioral data need to be put in there (data/eeg_data/ & data/behav_data)

The preprocessed EEG data and behavioral data can be downloaded here:
https://zenodo.org/records/13909754
The EEG data ("pp_reorder.mat") is the analysis-ready preprocessed version, namely the output of script 5 (s5_correctdata.m).
The behavioral data ("behav_res_pp.mat") is the output of script 1 (s0_extractdata.m).

Notes:
- The statistical values and results will not exactly match the reported results because there is randomness in the scripts (e.g., in classifier folding and label shuffling).
- Participant 16 did not finish the experiment and was excluded from all analyses.
- Two participants (2, 23) had poor data quality, and were excluded from EEG analyses, and one (7) was removed for trigger issues.
- At some point we switched EEG caps. The relevant changes are incorporated in the scripts.
- Make sure to change the path to your folders across scripts.
- If you want to check out the Psychopy scripts, ensure directories are lined up and you may need to use our version (2021.2.3; tip: you can run earlier versions from later versions)

Dependencies:
- FieldTrip
- MATLAB (R2020a used by us)
- Signal Analysis Toolbox
- MVPA Light
- Various scripts under /dependencies
- There may be some further dependencies not listed or not supplied under /dependencies for licensing reasons. If so, download the relevant files from the internet.
If you are unsure what package a function is from, search for "function_name.m matlab"

For questions, email s.....vanbree@gmail.com (first name)
