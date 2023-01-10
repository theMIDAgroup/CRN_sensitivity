# Statistical Sensitivity Index for CRN

Matlab code for performing a sensitivity analysis of a chemical reaction newtork (CRN) modeling cell signaling networks in colorectal cancer.

The performed analysis and the employed Statistical Sensitivity Index are described in the following manuscript.

G. Biddau, G. Caviglia, M. Piana and S.Sommariva (submitted) SSI: A Statistical Sensitivity Index for Chemical Reaction Networks in cancer


# For reproducing the results of the paper:

Code is written in Matlab R2017b.

Run main_correctness_index.m for testing the proposed SSI on the physiological network. 

Run main_Mutation_sensitivity_\*.m for testing the proposed SSI on a mutated network.


# External packages:

Our code makes use of the following tools:
* Robert (2021). symlog (https://github.com/raaperrotta/symlog), GitHub. Retrieved January 15, 2021.
* https://github.com/theMIDAgroup/CRC_CRN
