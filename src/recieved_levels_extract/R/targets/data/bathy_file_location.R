register_target(
  
  # Computer Directory Location for:
  # bathymetry file that will be used downstream
  # in the rejection sampler
  
  tar_target(
    
    name = bathy_file_location,
    command = {
      
      file.path(here::here('data/strm30_bathy.rds'))
    
      },
    format = "file")
             

)