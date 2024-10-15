register_target(
  
  # Here we do the sampling of the TOBL values
  # I need to read in the output from 
  # the RL manifest
  # the TOBL grids
  # One key difference is that the TOBL grids only need to be sampled 
  # once per cee
  
  tar_target(
    
    name = sample_ambient,
    command = {
      
      
      # Container to hold output
      rl_tobls <- vector(mode = 'list', length = length(calc_rl))
      # loop over one animal/cee combination at a time
      for(i in 1:length(calc_rl)){
        
        my_cee_id <- calc_rl[[i]][["cee_id"]]
        # Don't sample ambient for control animals
        if(my_cee_id %in% c("22_01", "22_02"))next()
        animal <- calc_rl[[i]][["pts_rl"]] %>% 
          sf::st_as_sf() 
          
        wind_path <- which(my_cee_id == process_wind_files[, "cee_id"])
        wind_file <- raster::raster(process_wind_files$wrast_path[wind_path])
        
        # get the TOBL Info
        animal$tobl <- raster::extract(wind_file, animal)
        
        # Write Animals to Disk - animal is both an sf object and data frame
        saveRDS(animal, file = here::here('data/proc_data', 
                                          paste0(unique(animal$deployid), "_", 
                                                 my_cee_id, '_tobl.rds')))
        
        # Store items in List 
        rl_tobls[[i]] <- list(deployid = unique(animal$deployid),
                            cee_id = my_cee_id,
                            pts_tobl = animal)
        
        
      } # end i loop over the list of each animal/cee combination
      
      rl_tobls|>
        purrr::compact()
      
      })
             

)