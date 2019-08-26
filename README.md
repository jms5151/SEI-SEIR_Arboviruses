# SEI-SEIR_Arboviruses
Weather-based mechanistic model of arbovirus transmission in Ecuador and Kenya.

# SEI-SEIR models:
Temperature-dependent model: SEI-SEIR_model_T.R
Temperature- and rainfall-dependent model: SEI-SEIR_model_TR.R
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_model_TRH.R
Temperature-, rainfall-, and humidity-dependent model with variation in c, T0, Tmax of temperature-dependent traits: 
        SEI-SEIR_model_with_trait_variation.R

# Data to initiate SEI-SEIR models:
Initial conditions (proportion of population in each compartment): LHS_inputs.csv
Posterior distribution of c, T0, Tmax for temperature-dependent traits: Random_sample_of_posterior_traits.csv

# Code to run model simulations:
Temperature-dependent model: SEI-SEIR_simulate_with_T_model.R
Temperature- and rainfall-dependent model: SEI-SEIR_simulate_with_TR_model_test_diff_rain_functions.R
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_simulate_with_TRH_model.R
Temperature-, rainfall-, and humidity-dependent model with variation in c, T0, Tmax of temperature-dependent traits: 
        SEI-SEIR_simulate_with_trait_variation.R

