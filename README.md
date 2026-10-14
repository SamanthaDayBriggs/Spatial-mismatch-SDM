GENERAL INFORMATION
This README.txt file was updated on January 30, 2026

Code and data for “Spatial mismatches of plant-pollinator interactions under future climate conditions”

Abstract: Global change could disrupt critical ecological interactions, especially if interacting species shift their geographic ranges at different rates or in different directions. Functional redundancy within generalized networks could make them more resilient to climate change than networks with highly specialized interactions. Here, we examined plant-pollinator interactions in the Rocky Mountain Iris (Iris missouriensis), which is distributed across broad climatic gradients over elevation and visited by a diversity of pollinators. Our ecological network analysis revealed that pollinator species richness and diversity is highest at mid-elevation sites. We then modeled range shifts in I. missouriensis and its pollinators, to test if climate change could generate spatial mismatches owing to species-specific changes in geographic distributions, disrupting plant-pollinator interactions. Our species distribution models of I. missouriensis and eleven of its pollinators revealed that future climate scenarios could reduce the probability of I. missouriensis co-occurrence with current pollinators. Thus, species-specific geographic responses to climate change could reshape the pollinator community even of charismatic plants with diverse pollinator assemblages.

Dates of data collection
2021-2022

ACCESS INFORMATION
Data and code will be made public upon publication. 

DATA & CODE FILE OVERVIEW
This data repository consist of 3 data files, 2 code scripts, and this README document, with the following data and code filenames and variables
Data files and variables
    1. poll_obs_year.csv [Column A: Elevation of data collection, B: Year of data collection, C-BH: Observed counts for each species]
    2. poll_obs_all.csv [Column A: Elevation of data collection, B-BG: Observed counts for each species]
    3. Iris_missouriensis.csv pulled directly from GBIF.org [Column A: datasetKey- ID of the GBIF dataset, B: basisOfRecord- nature of the recorded observation, C: occurrenceStatus- statement of the presence or absence of the species, D: species- species name, E: decimalLatitude- latitude of observation, F: decimalLongitude- longitude of observation, G: continent- continent of observation, H: stateProvince- state or province of observation, I: year- year of observation, J: month- month of observation, K: day- day of observation, L: issues- any interpretation issues, M: references- references associated with the occurrence, N: license- document for permission of use, O: geodeticDatum- the SRS of which the geographic coordinates given are based, P: gbifID- a unique ID for the occurrence in the GBIF database, Q: catalogNumber- an identifier for the record within the collection, R: elevation- elevation of observation, S: locality- a description of the place the occurrence was identified, T: county- the full name of the small administrative region that the occurrence was identified, U: type- nature of the resource, V: elevationAccuracy- value of error in elevation estimate, W: preparations- a list of any preparations or preservation methods used on the specimen]
Code scripts and workflow
    1. polldiversity.R, Rarefies yearly pollinator observation data then calculates Shannon Diversity and species richness. Uses GLMs to analyze patterns across elevations. Rarefies all data to calculate a beta diversity matrix to estimate species turnover across elevational sites. Generates the pollinator network across elevations. 
    2. SDM_overlap_elevation code.R, Provides a species distribution model workflow for Iris missouriensis and describes changes that should be made to repeat methods for each species. This pulls climate data directly from worldclim to generate the ensemble and binary models. GBIF occurrence files should be downloaded prior to running the code, but the Iris missouriensis occurrences are provided. The next step in the code overlays the binary models and calculates spatial overlaps under the different time and climate scenarios. The code then estimates the elevational range for each species based on a provided spatial extent. 

SOFTWARE VERSIONS
R v. 4.4.1
R packages:
vegan v. 2.6-8
car v. 3.1-3
stats v. 4.4.1
dplyr v. 1.1.4
tidyr v. 1.3.1
visreg v. 2.7.0
bipartite v. 2.23
raster v. 3.6-26
sp v. 2.1-4
megaSDM v. 2.0.0
sdm v. 1.2-55
geodata v. 0.6-2
ENMTools v. 1.1.2
dismo v. 1.3-14
earth v. 5.3.4
terra v. 1.7-78
gbm v. 2.2.2
kernlab v. 0.9-33
rJava v. 1.0-11
parallel v. 4.4.1
