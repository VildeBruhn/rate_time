# EFFECTS OF MODEL ADEQUACY AND SAMPLING ERROR ON EVOLUTIONARY RATES

__Article:__ 

__Authors:__ Vilde Bruhn Kinneberg<sup>*1</sup> and Kjetil Lysne Voje<sup>1</sup>

__Affiliation:__ <sup>1</sup>Natural History Museum, Univeristy of Oslo

__Contact:__ <sup>*</sup>v.b.kinneberg@nhm.uio.no

__Journal:__ 

__Year:__ 2024  

__Abstract:__ Evolutionary rates are found to decrease with increased time. How to link the high evolvability in microevolutionary processes to the slow evolution observed across macroevolutionary time is a conceptual challenge. In this study, we investigate the effects of model misspecification, sampling error and model identifiability on time scaling of evolutionary rates. We focus on the unbiased random walk and use both simulations and empirical data from the Phenotypic Time Series database (PETS). The rate estimated from the simulated unbiased random walk time series are shown to not scale with time, even when the time series are made incomplete and biased. When accounting for sampling error, model identifiability and model adequacy to avoid model misspecification, the negative scaling of rates with time still persists in the empirical time series. The results suggest that we have problems explaining the underlying trait dynamics in many empirical time series with current methods. The unbiased random walk model should be used with caution when comparing evolutionary rates across different time intervals, both in phenotypic time series and phylogenetic comparative analyses.

__Info:__ This repository contains scripts and data used for analyses in the publication.

__Responsibility:__ VBK is responsible for writing code and analysis of data. KLV has contributed with ideas and comments. The time series data are downloaded from the [PETS database](https://pets.nhm.uio.no/). 

__Files:__ 

_simulated_data –_ this folder conatins scripts and data used in the simulation section of the article. The main script is called simulations.R. The script uses the rest of the scripts and data in the folder. simulations.R is commented so that it should be possible to follow the instructions in the script to produce the results from the article.

_empirical_data –_ this folder conatins scripts and data used in the empirical section of the article. The two main scripts are called empirical.R and model_ident.R. These scripts use the rest of the scripts and data in the folder. empirical.R and model_ident.R are commented so that it should be possible to follow the instructions in the scripts to produce the results from the article. 

_supplementary_material_1 –_ The folder contains likelihood surface plots reffered to as SM1 in the article. The plots are produced in model_ident.R. The results are stored as zip-files.
