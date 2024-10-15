register_target(
  
  # Read in the cee data, cull out only the animals
  # that are the animals exposed to the scaled source
  # do some date formatting
  
  tar_target(
    name = cee_2019_depths,
    command = {

      cee_2019_depths <- readr::read_csv(cee_2019_depths_file_location, show_col_types = FALSE) %>% 
        mutate(pred_time = lubridate::ymd_hms(pred_time, tz = 'UTC')) %>% 
        filter(category != "depth_obs") %>% 
        mutate(cee_id = if_else(deployid == 'ZcTag093' & category == 'cee_first', '19_03', cee_id))
      
      return(cee_2019_depths)
    }
  )
  
)