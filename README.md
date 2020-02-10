# Software
All codes were created and run using R statistical software (https://www.r-project.org/) version 3.5.3 and Rstudio version 1.1.383 (https://rstudio.com/).

# SEI-SEIR_Arboviruses
This repository provides code to run climate-based mechanistic models of arboviral transmission. The models support of the manuscript "Climate explains geographic and temporal variation in mosquito-borne disease dynamics" by Caldwell et al (in prep). 

# SEI-SEIR models:
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_model_THR.R <br />
Intervention 1 - reduce contact rate between mosquitoes and people: SEI-SEIR_model_intervention_reduce_biteRate.R
Intervention 2 - reduce immature mosquito habitat: SEI-SEIR_model_intervention_reduceK.R
Intervention 3 - reduce mosquito abundance: SEI-SEIR_model_intervention_spray.R

# Data to initiate SEI-SEIR models:
Initial conditions (proportion of population in each compartment): LHS_inputs.csv <br />
Posterior distribution of c, T0, Tmax for temperature-dependent traits: Random_sample_of_posterior_traits.csv 
