register_target(
  
  # Locations for:
  # cee file that has the water column depths
  # for all the animals in 2019
  
  tar_target(
    
    name = cee_2019_depths_file_location,
    command = {
      
      here::here("data/2023-03-10_crawl-prediction-times.csv")
    
      },
    format = "file")
             

)