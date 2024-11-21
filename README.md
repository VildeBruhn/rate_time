# EFFECTS OF MODEL ADEQUACY, SAMPLING ERROR AND MODEL IDENTIFIABILITY ON EVOLUTIONARY RATES

__Article:__ Unpublished

__Authors:__ Vilde Bruhn Kinneberg<sup>1*</sup> and Kjetil Lysne Voje<sup>1</sup>

__Affiliation:__ <sup>1</sup>Natural History Museum, University of Oslo

__Contact:__ <sup>*</sup>v.b.kinneberg@nhm.uio.no

__Journal:__ NA

__Year:__ 2024  

__Abstract:__ It is challenging to compare evolutionary rates across clades because they depend on the time interval they are estimated over. In this study, we investigate the effects of model misspecification, sampling error, and model identifiability on the relationship between rates of phenotypic change and time in evolutionary time series. We first estimate rates of evolution as a parameter in the unbiased random walk model and use simulations to show that a lack of rate-time scaling is expected when the simulated data match the model used for rate estimation, even when time series are incomplete and biased. We then analyze 643 empirical evolutionary time series to assess whether accounting for model misspecification, sampling error, and model identifiability can reduce the negative scaling of rates with time. None of these measures appear to have an impact. We find a strong and consistent negative rate-time scaling even when removing rates estimated by models that inadequately describe trait dynamics in the data. Our results suggest that common models used in phylogenetic comparative studies and phenotypic time series analyses fail to accurately summarize trait evolution in empirical time series. This makes it difficult to compare estimated rates of evolution between datasets covering different time intervals.


__Info:__ This repository contains scripts and data used for analyses in the publication.

__Responsibility:__ VBK is responsible for writing code and analyses of data. KLV has contributed with ideas and comments. The time series data are downloaded from the [PETS database](https://pets.nhm.uio.no/). 

__Files__ 

_simulated_data –_ this folder conatins scripts and data used in the simulation section of the article. The main script is called simulations.R. The script uses the rest of the scripts and data in the folder. simulations.R is commented so that it should be possible to follow the instructions in the script to produce the results from the article.

_empirical_data –_ this folder conatins scripts and data used in the empirical section of the article. The five main scripts are called empirical.R, model_ident.R, empirical_>1Myr.R, empirical_segmented_reg.R and adeq_all.R. These scripts use the rest of the scripts and data in the folder. empirical.R, model_ident.R, empirical_>1Myr.R, empirical_segmented_reg.R and adeq_all.R are commented so that it should be possible to follow the instructions in the scripts to produce the results from the article. 

_supplementary_material_2 –_ The folder contains loglikelihood surface plots reffered to in the article. The plots are produced in model_ident.R. The results are stored as zip-files.
