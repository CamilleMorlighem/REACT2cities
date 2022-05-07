# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors 



This repository contains data and material to model and predict malaria risk (measured as *Pf*PR<sub>2-10</sub>) in two sub-Saharan African cities, Kampala (Uganda) and Dar es Salaam (Tanzania), using a set of environmental and socio-economic predictors derived from remote sensing imagery. These predictors are multi-resolution variables depicting the urban climate, the land use and the land cover. *Pf*PR<sub>2-10</sub> modelling and prediction are achieved using a popular machine learning algorithm, namely random forest (RF). 

<!-- Fore more information on the data and material presented here see this [paper](#to_be_added) (coming soon).  -->


## Data 

### Malaria prevalence data 

Malaria prevalence data (see folder `Data\Malaria_data\`) come from an open online malaria database recording survey data from several sources [[1](#1)], such as scientific papers or health surveys. The prevalence is measured as the *Plasmodium falciparum* Parasite Rate (i.e. the proportion of people infected by *Pf*) standardised over the two-to-ten age range (*Pf*PR<sub>2-10</sub>) [[2](#2)]. 

### Geospatial datasets

The following four geospatial datasets (see folder `Data\Variables\`) are used as covariates for modelling and mapping *Pf*PR<sub>2-10</sub> : 

* Pseudo-climate variables (COSMO) : Pseudo-climate variables are 1 km resolution raster grids produced by the COSMO-CLM model. They represent different climate variables as aggregate values (average, maximum or minimum) for the dry season (June to September 2014) [[3](#3)].
* Local climate zones (LCZ) : LCZ maps are binary maps of 100 m resolution classifying pixels into areas of uniform surface cover and structure with a specific temperature regime [[4](#4)]. They were derived applying a random forest (RF) classification algorithm to Landsat, USGS and Sentinel imagery from 2017 to 2019 [[5](#5)-[6](#6)]. 
* Land use (LU) and land cover (LC) : LC maps (0.5 m resolution) were derived from Pleiades satellite images acquired in 2013 for Kampala and in 2016 and 2018 for Dar es Salaam. The LC classification was performed using Computer Assisted Photo Interpretation, GEOBIA and machine learning [[7](#7)]. The LU maps (20 m resolution) were produced based on the LC maps and linear information extracted from OpenStreetMap [[8](#8)].
* Ancillary variables (Base): These are (i) the averaged NDVI and NDWI over period 2005-2019 derived from Landsat 5 and 8 satellite imagery (100 m resolution) and (ii) a SRTM digital elevation model of 2000 (30 m resolution) [[6](#6)]. 

### Other data

Besides malaria prevalence data and geospatial datasets, mapping *Pf*PR<sub>2-10</sub> requires a <a id="predgrid"></a>prediction grid of 1 km resolution (see folder `Data\Prediction_grid\`) covering the city extents and maps of the administrative boundaries (admin levels 4 and 5) of the cities (see folder `Data\Admin\`). 

## Material 

### Installing dependencies 

To run this code, you need a recent version of R (R 4.1+). Required R packages can be installed in R or R Studio with: 

```r
install.packages(c("spex", "sf", "sp", "raster", "stars", "XLConnect", "foreach", "dplyr", "plyr", "exactextractr", "doParallel", "iterators", "units", "corrplot", "ggplot2", "ggpubr", "mlr", "parallelMap", "pdp"))
```
Note that this code was developped using Windows and was not tested with macOS.

### Running the code

The whole process of modelling and mapping *Pf*PR<sub>2-10</sub> in Kampala and Dar es Salaam can be implemented by running the main script `Scripts\main_MRM.R`, which performs for each city the following steps : 

#### 1. Select malaria prevalence surveys

Through function `select.malaria.data`, the code selects surveys in the [database](#malaria-prevalence-data) (1) that were conducted over the 2005-2016 period, (2) that only include participants younger than 16 years old and (3) that are not Demographic and Health Surveys (DHS). You can change the following parameters for selecting malaria surveys : 

* `max_age` : integer defining maximum age of survey participants
* `survey_date` : vector with two integers respectively defining the start and end years of the considered time period
* `select_DHS` : boolean defining whether to select DHS or not 


#### 2. Extract covariates 

Using function `covariates.extraction`, covariates are extracted for the different input [geospatial datasets](#geospatial-datasets) either directly at the survey coordinates for 1 km resolution variables or in 1 km buffers for finer resolution variables. In addition, for a latter prediction stage, covariates are extracted for all pixel centroids of the 1 km resolution [prediction grid](#predgrid) covering the city extent. With the following parameter, you can choose whether to implement covariates extraction for the training and/or prediction stages: 

* `cov_extract_aims` : vector with `"training"` and/or `"prediction"` as elements

Optionally, you can use function `check.correlation` to create correlation matrices and check the correlation between all covariates and *Pf*PR<sub>2-10</sub>. 

#### 3. Create random forest (RF) models

*Pf*PR<sub>2-10</sub> is modelled using RF modelling via the function `rf.modelling`. RF models are built using a 10-repeated 5 folds spatial cross-validation (SCV) with 2 sub-folds for tuning the hyperparameters. Depending on the parameter `cov_selection_mode`, RF models are built using either all input covariates (`ALL`) or implementing a variable selection, i.e. a Recursive Feature Elimination (`RFE`). You can choose to implement RF modelling for two different purposes (defined by the parameter `variables_group`) : 
1) Build RF models with all input covariates altogether, i.e. Base, LULC, LCZ and COSMO (`variables_group = "all3Geo"`) ; 
2) Compare the importance of the different [geospatial datasets](#geospatial-datasets) for modelling *Pf*PR<sub>2-10</sub>, i.e. build RF models with Base, LULC, LCZ and COSMO separately (`variables_group = "GeoSpDt"`). 

In brief, you can change the following parameters for RF modelling: 

* `reps`, `folds`, `iters` : integers defining the number of repetitions, folds and sub-folds in the SCV 
* `cov_selection_mode` : character defining whether to implement a variable selection or not (`"ALL"` or `"RFE"`) 
* `variables_group` : character defining the purpose of RF modelling (`"all3Geo"` or `"GeoSpDt"`)


#### 4. Predict *Pf*PR<sub>2-10</sub>
 
Function `rf.prediction` is used to predict *Pf*PR<sub>2-10</sub> at 1 km resolution using the [prediction grid](#predgrid). The final predictive map is the average of the predictions made by the 50 models built in SCV (10 repetitions * 5 folds). Additional maps are created by aggregating the *Pf*PR<sub>2-10</sub> predictions by administrative boundaries. You can change the following parameter to define which RF models to use for the predictions: 

* `pred_mode` : character defining whether to use RF models implemented with a RFE (`"RFE-<i>"` where `i` is the i<sup>th</sup> RFE iteration) or not (`"ALL"`) 


## Results

Running this code creates a `Results` folder with four subfolders : (1) `Malaria_data`, (2) `Covariates`, (3) `RF_modelling` and (4) `Prediction`, which respectively contain the results of the [four steps previously described](#running-the-code). Interesting outputs are : 

* `Malaria_data\<city>\malaria_data_<city>.shp` : a shapefile storing the selected malaria surveys and their attributes 
* `Covariates\<city>\<city>_full_cov.csv` : a table storing all the extracted covariates used for training the RF models 
* `Covariates\<city>\Correlation\`: correlation matrices between covariates and the dependent variable 
* `RF_modelling\<variables_group>\<city>\<cov_selection_mode>\`: RF models, variables importance and dependency plots 
* `Prediction\<city>\final_plots\`: predictive maps at different resolution (1 km resolution and aggregated by administrative boundaries)

Here are examples of predictive maps of Dar es Salaam produced with these data and material. The left image shows the 1 km resolution predictions with the malaria survey data (black dots) and the right image shows the predictions aggregated by administrative areas (admin level 5). 

![alt_text](Pred_map.jpg)


## Credits

When using data from this repository, please cite the following data sources : 

| Dataset        | Source           
| :------------- |:-------------
| Malaria prevalence data       | [Snow, 2017a](#1)
| Pseudo-climate variables      | [Brousse *et al.*, 2020a](#3)      
| Local climate zones           |  [Brousse *et al.*, 2019](#5) (Kampala)<br>[Brousse *et al.*, 2020b](#6) (Dar es Salaam) |
| Land cover variables      | [Georganos and Grippa, 2020a](#9) (Kampala) <br> [Georganos and Grippa, 2020b](#10) (Dar es Salaam)   
| Land use variables      | [Georganos, 2020c](#11)
| Ancillary variables       | [Brousse *et al.*, 2020b](#6)

When using the methodology and/or code presented here, please cite [Morlighem, Chaiban *et al.*, 2022](#12).

## References

<a id="1"></a>[1] Snow, RW. The prevalence of Plasmodium falciparum in sub Saharan Africa since 1900. Harvard Dataverse. 2017. https://doi/10.7910/DVN/Z29FR0.

<a id="2"></a>[2] Snow, RW, Sartorius, B, Kyalo, D, Maina, J, Amratia, P, Mundia, CW, Bejon, P, & Noor, AM. The prevalence of Plasmodium falciparum in sub-Saharan Africa since 1900. Nature. 2017;550(7677):515-518. https://doi.org/10.1038/nature24059. 

<a id="3"></a>[3] Brousse, O, Wouters, H, Demuzere, M, Thiery, W, Van de Walle, J, & Lipzig, NV. The local climate impact of an African city during clear‐sky conditions—Implications of the recent urbanization in Kampala (Uganda). International Journal of Climatology. 2020;40(10):4586-4608. https://doi.org/10.1002/joc.6477. 

<a id="4"></a>[4] Stewart ID, Oke TR. Local Climate Zones for Urban Temperature Studies. Bulletin of the American Meteorological Society. 2012;93:1879–900. https://doi.org/10.1175/BAMS-D-11-00019.1. 

<a id="5"></a>[5] Brousse, O, Georganos, S, Demuzere, M, Vanhuysse, S, Wouters, H, Wolff, E, Linard, C and van Lipzig, NP. Urban climate using local climate zones in sub-Saharan Africa to tackle urban health issues. Urban Climate. 2019;27:227–242. https://doi.org/10.1016/j.uclim.2018.12.004.

<a id="6"></a>[6] Brousse, O, Georganos, S, Demuzere, M, Dujardin, S, Lennert, M, Linard, C, Snow, RW, Thiery, W, & van Lipzig, NPM. Can we use local climate zones for predicting malaria prevalence across sub-Saharan African cities? Environmental Research Letters. 2020;15:124051.  https://doi.org/10.1088/1748-9326/abc996. 

<a id="7"></a>[7] Grippa, T, Lennert, M, Beaumont, B, Vanhuysse, S, Stephenne, N, & Wolff, E. An Open-Source Semi-Automated Processing Chain for Urban Object-Based Classification. Remote Sensing. 2017;9(4):358. https://doi.org/10.3390/rs9040358. 

<a id="8"></a>[8] Grippa, T, Georganos, S, Zarougui, S, Bognounou, P, Diboulo, E, Forget, Y, Lennert, M, Vanhuysse, S, Mboga, N, & Wolff, E. Mapping Urban Land Use at Street Block Level Using OpenStreetMap, Remote Sensing Data, and Spatial Metrics. ISPRS International Journal of Geo-Information. 2018;7(7):246. https://doi.org/10.3390/ijgi7070246. 

<a id="9"></a>[9] Georganos, S, & Grippa, T. Kampala Very-High-Resolution Land Cover Map. Zenodo. 2020. https://doi.org/10.5281/zenodo.3711905.

<a id="10"></a>[10] Georganos, S, & Grippa, T. Dar Es Salaam Very-High-Resolution Land Cover Map. Zenodo. 2020. https://doi.org/10.5281/zenodo.3711903.

<a id="11"></a>[11] Georganos, S. Malaria in High-Resolution: Modelling and Mapping Plasmodium falciparum Parasite Rate using Very-High-Resolution Satellite Derived Indicators in Sub-Saharan African Cities. Zenodo. 2020. https://doi.org/10.5281/zenodo.3871497.

<a id="12"></a>[12] Morlighem, C*, Chaiban, C*, Georganos, S, Brousse, O, Vandewalle, J, van Lipzig, NPM, Wolff, E, Dujardin, S, Linard, C. Multi-satellite environmental and socio-economic predictors of vector-borne diseases in African cities: malaria as an example. 2022. * authors contributed equally.


