register_target(
  
  # Read in the cee data, 
  
  tar_target(
    name = cee_dat_hanks,
    command = {
      
      # Goniometer Data
      cee_dat_hanks <- readr::read_csv(cee_file_location, show_col_types = FALSE) %>% 
        dplyr::mutate(deployid = stringr::str_remove(deployid, "_DUML")) %>% 
        dplyr::filter(deployid %in% exposed_ans) %>% 
        dplyr::mutate(start_datetime = lubridate::ymd_hms(paste0(date_YYYYMMDD, start_HHMMSS))) %>% 
        dplyr::mutate(end_datetime = lubridate::ymd_hms(paste0(date_YYYYMMDD, end_HHMMSS)))
      
      
      return(cee_dat_hanks)
    }
  )
  
)