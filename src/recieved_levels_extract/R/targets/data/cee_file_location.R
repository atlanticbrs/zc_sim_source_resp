register_target(
  
  # Locations for:
  # cee file that has the animal by cee information
  
  tar_target(
    
    name = cee_file_location,
    command = {
      
      file.path("data/tagstreams_by_cee_details_filt.csv")
    
      },
    format = "file")
             

)