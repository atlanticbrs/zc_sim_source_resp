register_target(
  
  # Read in the cee data, cull out only the animals
  # that are the animals exposed to the scaled source
  # do some date formatting
  
  tar_target(
    name = cee_dat,
    command = {
      
      # Goniometer Data
      cee_dat <- readr::read_csv(cee_file_location, show_col_types = FALSE) %>% 
        dplyr::mutate(deployid = stringr::str_remove(deployid, "_DUML")) %>% 
        dplyr::filter(deployid %in% exposed_ans,
                      cee_type %in% c("simulated mfas", "control")) %>% 
        dplyr::mutate(start_datetime = lubridate::ymd_hms(paste0(date_YYYYMMDD, start_HHMMSS))) %>% 
        dplyr::mutate(end_datetime = lubridate::ymd_hms(paste0(date_YYYYMMDD, end_HHMMSS)))
      
      
      return(cee_dat)
    }
  )
  
)