register_target(
  
  # Read in the RL Manifest data, that has the RLs for all
  # the scaled source animals
  # you will also need this when calculating distance to ship
  # and when calculating RLs
  
  tar_target(
    name = rl_dat,
    command = {
      
      rl_dat <- readr::read_csv(rl_times_file_location,
                                show_col_types = FALSE)
      
      return(rl_dat)
    }
  )
  
)