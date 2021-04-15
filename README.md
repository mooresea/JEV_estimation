# JEV_estimation

This repository contains data and code to recreate the analysis for the MedRxiv preprint: "The current burden of Japanese encephalitis and the estimated impacts of vaccination: Combining estimates of the spatial distribution and transmission intensity of a zoonotic pathogen."

https://medrxiv.org/cgi/content/short/2021.04.08.21255086v1

The script <estimate_foi_stan_age_groups.R> will estimate annual FOI and past vaccination coverage estimates for each country from age-structured incidence data

The script <vimc_calc_impact_present_2020.R> uses the FOI estimates, and baseline estimates of the number of people at risk of infection in each country, and vaccination coverage estimates to estimate the annual number of cases and deaths for each country under current vaccination coverage levels and a "no vaccination" counterfactual scenario

The scripts <vimc_calc_impact_present_2020_low.R> and <vimc_calc_impact_present_2020_hi.R> modify the above analysis using low and high population-at-risk estimates
