# zc_sim_source_resp

This respository contains the data and codue used to produce analyses and figures for the accepted version of:

> Southall BL, Schick RS, Cioffi WR, DeRuiter SL, Foley HJ, Harris CM, Harshbarger AE, Joseph JE, Margolina T, Nowacek DP, Quick NJ, Swaim ZT, Thomas L, Waples DM, Webster DL, Wisse JH, Read AJ. 2025. Behavioral responses of goose-beaked whales (*Ziphius cavirostris*) to simulated military sonar. Ecosphere. doi: https://doi.org/10.1002/ecs2.70501

If you have questions about anything in this repository or are interested in using this daaset please contact Brandon Southall (brandon.southall@sea-inc.net).

## details
- directories:
- `00_data_input` contains the metadata on timing of sound exposure treatments
- `src` contains the source code and is organized into directories for each analysis or figure. Input data which is specifically relevant to each of these analyses is included in these directories.
- `01_shared_data_products` contains intermediate data files created by code in `src` and shared across different analyses.
- To run some analyses you will need to download `02_large_data` (see below).


## included source code directories
- gam -- generalized linear models (dive response)
- horizontal displacement
- multiscale -- combining dtag and stag data
- pre post -- dive response
- pre post explore -- dive response data exploration
- received levels extract -- calculated received levels of exposures
- sde -- dtag dive analysis
- social -- group composition analysis before during and after exposures
- univariate -- distributions of dive patterns compared to exposure dives

## large input data files
Some of the directories require very large input data files to run that don't fit on github. They can be downloaded at this link: https://duke.box.com/v/zcss-large-data. The contents can be deposited into `02_large_data`. All other required input data files should already be included in the repository

## dependencies

Source code should run on R 4.4, but a couple of dependencies are not available from CRAN. See below for some installation possibilities.

### `rgdal` can be installed on windows for R 4.4 with the folowing
```r 
url <- "https://download.r-forge.r-project.org/bin/windows/contrib/4.4/rgdal_1.6-7.zip"
  install.packages(url, type="source", repos=NULL)
```

### `colorblindr` requires non-CRAN versions of `cowplot` and `colorspace`
```r
remotes::install_github("wilkelab/cowplot")
install.packages("colorspace", repos = "http://R-Forge.R-project.org")
remotes::install_github("clauswilke/colorblindr")
```

### `crawlUtils`
```r
install.packages('crawlUtils', repos=c('https://dsjohnson.r-universe.dev','https://cloud.r-project.org'))
```

## targets
recieved levels extract and horizontal displacement use the package `targets`. Start by sourcing the `_targets.R` file. `targets::tar_visnetwork()` will give you an idea of the workflow. `targets::tar_make()` will run the entire workflow (make sure the large data files are already downloaded and deposited in the correct directory (see above). See more details at: https://books.ropensci.org/targets/.

# citation

Please cite this dataset as:

> Southall BL, Schick RS, Cioffi WR, DeRuiter SL, Foley HJ, Harris CM, Harshbarger AE, Joseph JE, Margolina T, Nowacek DP, Quick NJ, Swaim ZT, Thomas L, Waples DM, Webster DL, Wisse JH, Read AJ. zc_sim_source_resp: Dataset and code: Behavioral responses of goose-beaked whales (*Ziphius cavirostris*) to simulated military sonar. doi: https://doi.org/10.5281/zenodo.17640059

Please cite the manuscipt as:

>  Southall BL, Schick RS, Cioffi WR, DeRuiter SL, Foley HJ, Harris CM, Harshbarger AE, Joseph JE, Margolina T, Nowacek DP, Quick NJ, Swaim ZT, Thomas L, Waples DM, Webster DL, Wisse JH, Read AJ. 2025. Behavioral responses of goose-beaked whales (*Ziphius cavirostris*) to simulated military sonar. Ecosphere. doi: https://doi.org/10.1002/ecs2.70501

# contact information

If you have questions about anything in this repository or are interested in using this dataset please contact Brandon Southall (brandon.southall@sea-inc.net).
