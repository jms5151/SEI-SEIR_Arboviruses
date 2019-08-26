# SEI-SEIR_Arboviruses
This repository provides code to run weather-based mechanistic models of arboviral transmission. The models support of the manuscript "Weather explains geographic and temporal variation in mosquito-borne disease dynamics" by Caldwell et al (in prep). 

# SEI-SEIR models:
Temperature-dependent model: SEI-SEIR_model_T.R <br />
Temperature- and rainfall-dependent model: SEI-SEIR_model_TR.R <br />
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_model_TRH.R <br />
Temperature-, rainfall-, and humidity-dependent model with variation in c, T0, Tmax of temperature-dependent traits: SEI-SEIR_model_with_trait_variation.R

# Data to initiate SEI-SEIR models:
Initial conditions (proportion of population in each compartment): LHS_inputs.csv <br />
Posterior distribution of c, T0, Tmax for temperature-dependent traits: Random_sample_of_posterior_traits.csv 

# Code to run model simulations:
Temperature-dependent model: SEI-SEIR_simulate_with_T_model.R <br />
Temperature- and rainfall-dependent model: SEI-SEIR_simulate_with_TR_model_test_diff_rain_functions.R <br />
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_simulate_with_TRH_model.R <br />
Temperature-, rainfall-, and humidity-dependent model with variation in c, T0, Tmax of temperature-dependent traits: SEI-SEIR_simulate_with_trait_variation.R

