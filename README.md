# Software
All codes were created and run using R statistical software version 3.5.3 and Rstudio version 1.1.383.

# About SEI-SEIR_Arboviruses
This repository provides code to run climate-based mechanistic models of arboviral transmission. The models support of the manuscript "Climate predicts geographic and temporal variation in mosquito-borne disease dynamics" by Caldwell et al. Nature Communications (accepted). 

# models
Temperature-, rainfall-, and humidity-dependent model: SEI-SEIR_model_THR.R <br />
Temperature-, rainfall-, and humidity-dependent model with trait variation: SEI-SEIR_model_with_trait_variation.R <br />
Code to set up simulations: SEI-SEIR_simulation_setup.R <br />

# data
Initial conditions (proportion of population in each compartment): LHS_inputs.csv <br />
Posterior distribution of c, T0, Tmax for temperature-dependent traits: Random_sample_of_posterior_traits.csv  <br />
Epidemic characteristics datasets: (1) Duration_intensity_outbreaks.csv and (2) Number_outbreaks.csv <br />
Time series of temperature, humidity, and rainfall for study sites: climate_data.RData <br />
Time series of vector data for study sites: vector_data.csv. Only adult Ae. aegypti data was collected at Ecuador sites. <br />
Time serires of arboviral data for study sites: arbovirus_data.csv. Data is only available for Kenya sites.  <br />

# simulations
A single code is provided to run each model with appropriate data.

# analyses 
Compare model output with observational data: pairwise_correlation.R <br />
Compare model output with different socio-ecological characteristics: socio_ecological_relationships.R <br />
Compare predictions and observations for different epidemic characteristics: epidemic_characteristics.R <br />
