# AVO_earthquake_classification
Code to accompany Power and Roman (2024) - calculation of various metrics and their use for volcanic earthquake autoclassification

Last updated August 18, 2024 by Diana Roman

Code to accompany Power and Roman (2024), JVGR. This code runs on pairs of mseed and pha files with a shared name (as in the three examples).
It selects the three closest stations to the event epicenter and uses these to calculate a stacked/average spectrogram around the earthquake (based on the arrival times) and then calculates Frequency Index, peak frequencies, kurtosis and skewness (quality critera) and then automatically assigns a class (Low-Frequency or High-Frequency) to the event.

Note: You may need to create MLtrue.mat table for accurate magnitude reporting in output files (see example). This is only if the .pha files does not contain the correct magnitude. You can ignore this if it is not needed. 

To run the code: 

1) First, set the appropriate variables (frequency bands for FI calculation, bandpass filter corners) in autoclass.m before running

2) Place .pha and .mseed pairs in the directory with the matlab code files and run as
> autoclass

3) When prompted select mode: 
1 - analysis
2 - manual checking (displays automatic classification result)
3 - manual checking (blind - only shows stacked spectrogram and top picks)

It will return an output in the form of the varible 'Table'

If run in mode 1, the final column (manual check) will be 0
In modes 2 or 3, the final column (manual check) is as follows:
1 - High-Frequency
2 - Low-Frequency
3 - Uncertain
