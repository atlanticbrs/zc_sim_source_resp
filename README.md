# For review only!

## included source code directories
- gam -- generalized linear models (dive response)
- horizontal displacement
- multiscale -- combining dtag and stag data
- pre post -- dive response
- pre post explore -- dive response data exploration
- recieved levels extract -- calculated recieved levels of exposures
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
