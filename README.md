# tDCS Motor Learning MEP TEP

This is the code used for our study on the dose response of high-intensity tDCS on motor learning, where we measured motor learning behavior and changes in corticospinal (MEP) and cortico-cortico (TEP) excitability with TMS. For more details, read our manuscript at: [https://doi.org/10.1162/imag_a_00431](https://doi.org/10.1162/imag_a_00431)

An archive of our data, manuscripts, and code for this project can be found on OSF: [https://osf.io/gyfae/](https://osf.io/gyfae/)

This code was last tested using MATLAB R2023a and may not be compatible with other versions. Version R2019b and later required.

The function `addeeglab` is shared across processing and analysis code.

## Processing

These are the pipelines used to process MEP/EMG and TEP/EEG data prior to analysis.

### Requirements

Make sure these are on the MATLAB path!

- **Libeep MATLAB importer** (last tested with 3.3.177)
- **EEGLAB** (last tested with 2024.0)
- **TESA plugin for EEGLAB** (last tested with 1.1.1)
- **BVA-iO plugin for EEGLAB** (last tested with 1.73)
- **FastICA** (last tested with 2.5)


### EMG pipeline

1. Manually enter trigger information into TMStrials.csv, noting the first and last trigger numbers for each hand and TMS session
2. Some recordings may contain extra/false TMS triggers. Check for extra triggers using `checkTMSTriggers`
3. Extract MEP data using `getMEPprepost`
4. Manually add newly processed data to `MEPdata.mat`
5. Visualize MEPs one subject at a time using `viewMEPprepost`

In case of trigger errors in `read_eep_trial` (Libeep), check trigger numbers using `read_eep_trg` (also from Libeep)

### EEG pipeline

1. Run the script `loopSaveEpochs` to detect and extract TMS epochs from raw EEG recordings, looping through multiple subjects. This outputs a table `TT` that shows how many epochs are detected. Check that the number is correct for each subject
2. If extra or missing trials are found, inspect recordings individually using `inspectTriggers`. This function also overlays detected TMS pulses with true TMS triggers recorded with EMG
3. Remove any extra triggers using `removeEpochs`
4. Review epoched data using `inspectEpochs`
5. Rerun `saveEpochs` and repeat inspection steps for individual subjects if necessary
6. Run the script `processEEG` to run FastICA and TESA artifact removal for all subjects
7. Run `plotTEP` to inspect averaged TEPs for individual subjects
8. Run the script `saveTEPs` to save all processed TEP data into one .mat file

## Analysis

At the end of processing, we have MEP and TEP data compiled into two .mat files: `MEPdata.mat` and `TEPdata.mat`. Behavioral data (motor learning performance and typing speed) are not processed. All of these data are combined into one file, `combinedData.mat`, for analysis.

Additional data files are required to run the analysis code, including data on group sorting, demographics, and sensation rating questionnaires. See the data on our OSF repository as an example.

### Analysis pipeline

1. Run the script `combineData` to combine behavioral, MEP, and TEP data into a single file.
2. Run the scripot `figsandstats` for final analysis and visualization. All figures and statistics in the manuscript are created from this script.

### Other scripts
- `TT_analysis` - used to sort subjects into groups based on typing speed. See manuscript for details.
- `SimulateData` - used to simulate an effect across groups, using a copy of existing data. We used this to check that our analysis code was able to correctly determine an effect or lack thereof.

## Task

This is the code used to administer our motor sequence learning task. It is run as an app, installed from `SequenceTask2_App.mlappinstall`. The app is run from the script `runFTT`, which prompts the user to enter the subject ID as well as task sesssion information.

## Acknowledgements

Davide Bonfanti, for providing reference code for the EEG processing pipeline

Zhenous Hadi Jafari, for writing the code to display sensation rating questionnaire results

## License
GNU General Public License (GPL) 3.0

See LICENCE.md for details.

Copyright (c) 2025 Gavin Hsu, Parra Lab, CCNY, New York, NY