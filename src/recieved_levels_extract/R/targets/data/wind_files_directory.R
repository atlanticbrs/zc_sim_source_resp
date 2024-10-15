register_target(
  
  # Locations for: Wind files for all CEEs
  # These will be used to store the right wind raster to sample 
  
  tar_target(
    
    name = wind_files_directory,
    command = {
      
      wind_files <- dir(here::here("../../02_large_data/wind"), 
                      pattern = "^Wind.*Noise.*nc$",
                      recursive = TRUE,
                      full.names = TRUE)
      
      wind_files
      
      })
             

)