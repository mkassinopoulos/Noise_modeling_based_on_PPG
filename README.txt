
Michalis Kassinopoulos, PhD
Graduate Program in Biological and Biomedical Engineering, McGill University, Montreal, Canada
mkassinopoulos@gmail.com

Date: 12-August-2021

==================================================================

This repository provides scripts used for our work presented in: 
Kassinopoulos, M., Mitsis, G.D., 2021. Physiological noise modeling in fMRI based on the pulsatile component of photoplethysmograph. Neuroimage 118467. https://doi.org/https://doi.org/10.1016/j.neuroimage.2021.118467


Abstract: The blood oxygenation level-dependent (BOLD) contrast mechanism allows someone to non-invasively probe changes in deoxyhemoglobin content. As such, it is commonly used in fMRI to study brain activity since levels of deoxyhemoglobin are indirectly related to local neuronal activity through neurovascular coupling. However, the BOLD signal is severely affected by physiological processes as well as motion. As such, several noise correction techniques have been developed through the years to correct for the associated confounds. This study sought to refine model-based techniques that utilize the photoplethysmograph (PPG) signal. RETROICOR, a technique commonly used to model fMRI fluctuations induced by cardiac pulsatility was compared with a modified version of it, named cardiac pulsatility model (CPM) that is based on convolution filtering. Further, this study investigated whether the variations in the amplitude of the PPG pulses (PPG-Amp) covary with variations in amplitude of pulse-related fMRI fluctuations as well as with systemic low frequency oscillations (SLFOs) present in the global signal (i.e. mean fMRI timeseries averaged across all voxels in gray matter). Capitalizing on 3T fMRI data from the Human Connectome Project, CPM was found to explain significantly more variance in fMRI compared to RETROICOR particularly for subjects that presented high variance in heart rate during the scan. The amplitude of the fMRI pulse-related fluctuations did not seem to covary with PPG-Amp. That said, PPG-Amp explained significant variance in the GS that did not seem to be attributed to variations in heart rate or breathing patterns. In conclusion, our results suggest that the techniques proposed here can model high-frequency fluctuations due to pulsation as well as low-frequency physiological fluctuations more accurately than model-based techniques typically employed in fMRI studies.  

====================================================================

The User guide will be available shortly and will explain what each script is related to and how to reproduce some of the results presented in the study.

Please do not hesitate to contact me if you have any questions related to the use of these scripts.

====================================================================
License

The scripts in this repository are licensed under the Apache License, Version 2.0. See LICENSE.md for the full license text.
