# Software
All codes were created and run using R statistical software version 3.5.3 and Rstudio version 1.1.383.

# About SEI-SEIR_Arboviruses
This repository provides code to run climate-based mechanistic models of arboviral transmission. The models support of the manuscript "Climate explains geographic and temporal variation in mosquito-borne disease dynamics" by Caldwell et al (in prep). 

# models
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_model_THR.R <br />
Temperature-, rainfall-, and humidity-dependent model with trait variation: SEI-SEIR_model_with_trait_variation.R <br />
Code to set up simulations: SEI-SEIR_simulation_setup.R <br />

# data
Initial conditions (proportion of population in each compartment): LHS_inputs.csv <br />
Posterior distribution of c, T0, Tmax for temperature-dependent traits: Random_sample_of_posterior_traits.csv 
Additional data will be upload upon peer-reviewed publication

# simulations
A single code is provided to run each model with appropriate data.

# analyses 
Compare model output with observational data: pairwise_correlation.R <br />
Compare model output with different socio-ecological characteristics: socio_ecological_relationships.R <br />
Compare predictions and observations for different epidemic characteristics: epidemic_characteristics.R <br />
