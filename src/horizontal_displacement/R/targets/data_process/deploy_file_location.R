register_target(
  
  # Locations for:
  # the deployment file that I use for the first location of each animal
  
  tar_target(deploy_file_location,
             here::here("../../01_shared_data_products/filter_sattag/metadata.csv"),
             format = "file")

)