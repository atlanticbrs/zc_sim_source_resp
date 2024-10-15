register_target(
  # Create the Unique Time Vector at the 5-minute intervals
  # for each animal/CEE combination
  tar_target(
    name = cee_times,
    command = {
      
      cee_dat_unique <- cee_dat %>%
        dplyr::group_by(cee_id) %>% 
        dplyr::slice_head()
      
      my_cee_times <- vector(mode = "list", length = nrow(cee_dat_unique))
      
      
      for(i in 1:nrow(cee_dat_unique)){
        
        # ID of cee
        my_id <- cee_dat_unique$cee_id[i]
        
        # Logic here is that we have different types of cees
        # and the controls have no sound, so we can't peel the
        # temporal information from the netCDF files
        if(cee_dat_unique$cee_type[i] == 'control'){
          
          cee_01_seq <- as.POSIXct(seq(cee_dat_unique$start_datetime[i], 
                                       cee_dat_unique$end_datetime[i], 
                                       by = "5 mins"), 
                                   origin = '1970-01-01', tz = 'UTC')

          if(max(cee_01_seq) != cee_dat_unique$end_datetime[i]){
            
            cee_01_seq <- c(cee_01_seq, cee_dat_unique$end_datetime[i])
            
          }
          
          rvec <- lubridate::round_date(cee_01_seq, '5 mins')
          rvec <- rvec[rvec > min(cee_01_seq) & rvec < max(cee_01_seq)]
          cee_tvec <- c(cee_01_seq[1], rvec, cee_01_seq[length(cee_01_seq)])
          
        } else {
          
          cee_tvec <- rl_dat %>% 
            filter(cee_id == my_id) %>% 
            pull(rl_time)          
          
        }

        
        my_cee_times[[i]] <- list(cee_id = my_id, 
                                  times = cee_tvec)
        
      } # end loop over cees to extract the time sequences
      
      names(my_cee_times) <- cee_dat_unique$cee_id
      
      return(my_cee_times)
    }
  )

)