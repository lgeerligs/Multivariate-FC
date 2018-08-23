# Multivariate-FC
Example scripts for computing functional connectivity using distance correlation.

Based on the paper:
Geerligs L., Cam-CAN & Henson R.N. (2016) Functional connectivity and structural covariance between regions of interest can be measured more accurately using multivariate distance correlation. Neuroimage. 2016; 135: 16-31
https://www.sciencedirect.com/science/article/pii/S1053811916300878?via%3Dihub 

The **dcor_uc** script will compute the U-centered (bias corrected) distance correlation between two input matrices. 
The **dcor_dc** script will compute the double centered (non bias corrected) distance correlation between two input matrices. 

The **example_connectivity_analysis** script is a demonstration of how to analyze functional connectivity data using distance correlation. 
This script loads example data that can be downloaded from http://imaging.mrc-cbu.cam.ac.uk/imaging/Geerligs_DistCor 

