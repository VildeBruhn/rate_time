# RATE-TIME SCALING IN PHENOTYPIC EVOLUTION: LIMITATIONS OF CURRENT MODELS IN CAPTURING TEMPORAL DYNAMICS

__Article:__ Unpublished

__Authors:__ Vilde Bruhn Kinneberg<sup>1*</sup> and Kjetil Lysne Voje<sup>1</sup>

__Affiliation:__ <sup>1</sup>Natural History Museum, University of Oslo

__Contact:__ <sup>*</sup>v.b.kinneberg@nhm.uio.no

__Journal:__ Evolution

__Year:__ 2025  

__Abstract:__ Evolutionary rates correlate negatively with time, which makes it complicated to compare rates across lineages that diversified on different time intervals. The causes of this correlation are debated. In this study, we use evolutionary time series to investigate the effects of model misspecification, sampling error, and model identifiability on the relationship between rates of phenotypic change and time. Using simulations, we first show that rates of evolution estimated as a parameter in the unbiased random walk (Brownian motion) model lack a rate-time scaling when data has been generated using this model, even when time series are made incomplete and biased. This indicates that it is theoretically possible to estimate rates that are not time correlated. We then analyze 643 empirical time series to assess whether accounting for model misspecification, sampling error, and model identifiability reduce the negative scaling, but none appear to have a significant impact. This suggests that the rate-time correlation requires an explanation grounded in evolutionary biology, and that common models used in phylogenetic comparative studies and phenotypic time series analyses often fail to accurately describe trait evolution in empirical data. Making meaningful comparisons of estimated rates between clades and lineages covering different time intervals remains a challenge.

__Info:__ This repository contains scripts and data used for analyses in the publication.

__Responsibility:__ VBK is responsible for writing code and analyses of data. KLV has contributed with ideas and comments. The time series data are downloaded from the [PETS database](https://pets.nhm.uio.no/). 

__Files__ 

_simulated_data –_ this folder conatins scripts and data used in the simulation section of the article. The main script is called simulations.R. The script uses the rest of the scripts and data in the folder. simulations.R is commented so that it should be possible to follow the instructions in the script to produce the results from the article.

_empirical_data –_ this folder conatins scripts and data used in the empirical section of the article. The two main scripts are called empirical.R and model_ident.R. These scripts use the rest of the scripts and data in the folder. empirical.R and model_ident.R are commented so that it should be possible to follow the instructions in the scripts to produce the results from the article. 

_supplementary_material_2 –_ The folder contains log-likelihood surface plots reffered to in the article. The plots are produced in model_ident.R. The results are stored as zip-files.
