# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(future)
# library(future.callr)
# library(crew)
#
# tar_option_set(
#   controller = crew::crew_controller_local(workers = 2)
# )


# Set target options:
tar_option_set(
  packages = c("tibble",
               "tidyverse",
               "data.table",
               "lubridate",
               "tidyr",
               "here",
               "sf",
               "crsuggest",
               "trip",
               "crawl",
               "crawlUtils",
               "ncdf4",
               "raster",
               "rgdal",
               "ctmcmove",
               "mgcViz",
               "mgcv",
               "colorblindr"), # packages that your targets need to run
  format = "rds", # default storage format
  seed = 314
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# project target container
target_list = list()

# add a tar_target object to a target container
#
# Parameters:
#   target - object to add to target_container
#   target_container - list() object with all targets
#   env - location in which target_container is defined
register_target = function(
    target, target_container = target_list, env = parent.frame()
) {
  # deparse target_container to a character, if needed
  if(is.character(target_container)) {
    tgt = target_container
  } else {
    tgt = deparse(substitute(target_container))
  }
  # append target to the container; update container
  assign(x = tgt, value = c(get(x = tgt), target), envir = env)
}

## load R files and workflows; must use register_target to populate target_list
lapply(list.files("R", full.names = TRUE, recursive = TRUE, pattern = '\\.R'), 
       source)

# pass on the list of updated targets to workflow manager
target_list
