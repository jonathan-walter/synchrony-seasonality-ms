# synchrony-seasonality-ms
Data and code accompanying Walter et al. (2023) Seasonality in environment and population processes alters spatial synchrony. The American Naturalist. doi: TBA. This study uses theoretical models to explore how seasonality influences population spatial synchrony.

There are four main code files and several helper code files that contain functions for simulation models.
The main code files are analytical_results.R, simulation_studies.R, sim_analytical_comparison.R, and cross_synchrony.R.
analytical_results.R explores our analytical solution and reproduces Figures 1 and 2.
simulation_studies.R runs simulation studies and reproduces Figures 3-5 and S2-S11.
sim_analytical_comparison.R examines consistency between our approximate analytical solution and simulation studies where some assumptions facilitating a tractable analytical solution are relaxed. This script reproduces Figure S1.
cross_sychrony.R assesses evidence for synchrony between environmental drivers (average temperature and total precipitation) in different seasons (winter and summer). It reproduces results reported in Online Appendix B.

The helper code files are simmod_main.R, simmod_main_deterministic.R, simmod_alt_nowinter.R, simmod_alt_sameenv.R, and simmod_main_arEnv.R.
simmod_main.R is the main formulation of the simulation model presented in the main text.
simmod_main_deterministic.R is a deterministic (no environmental noise) version of the main model that's used in some intermediate steps, e.g. determining whether population growth is undercompensatory or overcompensatory.
simmod_alt_nowinter.R is an alternate formulation of the main model used for comparison with a model with no overwintering phase.
simmod_alt_sameenv.R is an alternate formulation of the main model used for comparison with a model where both the breeding and overwintering phases experience the same environmental conditions.
simmod_main_arEnv.R is an elaboration of the main model where the environmental noises have autoregressive (AR) temporal structure. It is used to produce Figures S10 and S11.

There are four data files that facilitate the exploration of cross-season synchrony between environmental drivers.
These data are taken from the PRISM [https://prism.oregonstate.edu] interpolated (4 km grid) climate dataset at 1000 randomly selected geographic coordinates in the conterminous United States and reflect time series of seasonally aggregated (mean for temperature, sum for precipitation) temperature and precipitation from 1990 to 2009 (20 years).
Summer was June-July-August and Winter was December-January-February.
Each is in wide format with each row corresponding to a location and each column corresponding to a year.
