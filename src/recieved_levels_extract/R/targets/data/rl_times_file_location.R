register_target(
  
  # Locations for:
  # the deployment file that I use for the first location of each animal
  
  tar_target(rl_times_file_location,
             here::here("data/scaled_source_rl_manifest.csv"),
             format = "file")

)