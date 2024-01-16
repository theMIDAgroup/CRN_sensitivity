# Statistical Sensitivity Index for CRN

Matlab code for performing a sensitivity analysis of a chemical reaction newtork (CRN) modeling cell signaling networks in colorectal cancer.

The performed analysis and the employed Statistical Sensitivity Index are described in the following manuscript.

G. Biddau, G. Caviglia, M. Piana and S. Sommariva (submitted) SSI: A Statistical Sensitivity Index for Chemical Reaction Networks in cancer


# For reproducing the results of the paper:

Code is written in Matlab R2017b.

Run main_correctness_index.m for testing the proposed SSI on the physiological network. 

Run main_Mutation_sensitivity_\*.m for testing the proposed SSI on a mutated network.

* Run main_correctness_index.m for testing the proposed SSI on the physiological network end on the subnetwork related to p53.

* Run main_mutation_sensitivity.m for testing the proposed SSI on a mutated network. Four different networks were tested each one including a different mutation, namely a gain of function mutation (GoF) of gene KRAS, and a loss of function (LoF) of gene APC, SMAD4, and TP53.

* Run main_add_drugs_sensitivity.m for testing the proposed SSI on network affected by a GoF mutation of KRAS and considering the combined action of two drugs, namely Dabrafenib (DBF) and Trametinib (TMT).

All tests are based on data concerning a colorectal cells during their G1-S phase. The data of the network are collected in the folder ./data. Specifically, data are provided as MATLAb structures divided in four fields:
*	Species contains information about the chemical species, such as names and initial standard values
*	Rates contains information about the rate constants’ names and value
*	Reactions 
*	Matrix contains all the matrices and vectors described in the article, such as S, N and v.


# External packages:

Our code makes use of the following tools:
* Robert (2021). symlog (https://github.com/raaperrotta/symlog), GitHub. Retrieved January 15, 2021.
* https://github.com/theMIDAgroup/CRC_CRN
