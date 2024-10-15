register_target(
  
  # Here I make the depth data frame that holds the minimum/max
  # depth for each animal
  
  tar_target(
    
    name = an_series_df,
    command = {
      
      # Get the Depths correct
      series_cee_dat <- cee_dat %>% 
        filter(ser == TRUE,
               is.na(ser_gaps)) 
      
      my_depths <- vector(mode = 'list', length = nrow(series_cee_dat))
      for(i in 1:nrow(series_cee_dat)){
        
        animal <- series_cee_dat$deployid[i]
        cee_id <- series_cee_dat$cee_id[i]
        my_cee <- which(cee_id == names(cee_times))
        my_cee_times <- cee_times[[my_cee]]$times

        
        if(cee_id %in% c("19_01", "19_02", "19_03", "19_04")){
          
          this_an_series <- cee_2019_depths %>% 
            dplyr::filter(deployid == animal) %>% 
            dplyr::filter(pred_time %in% my_cee_times) %>% 
            dplyr::select(deployid, cee_id, datetime = pred_time, depth_min, depth_max) %>% 
            as.data.frame()
          
        } else {
          
          
          # need all of the data to do the forward/backward padding
          this_an <- series_data %>% 
            dplyr::filter(deployid == animal) %>% 
            as.data.frame()
          
          this_an_series <- series_data %>% 
            dplyr::filter(deployid == animal) %>% 
            dplyr::filter(datetime %in% my_cee_times) %>% 
            dplyr::mutate(depth_min = depth - d_range,
                          depth_max = depth + d_range,
                          cee_id = cee_id) %>% 
            dplyr::select(deployid, cee_id, datetime, depth_min, depth_max) %>% 
            as.data.frame()
          
          this_an_series$depth_min <- if_else(this_an_series$depth_min < 0, 0, this_an_series$depth_min)
          
          # Pad the first and last values
          # and interpolate their depths
          first_depth_time <- my_cee_times[1]
          last_depth_time <- my_cee_times[length(my_cee_times)]
          fidx <- max(which(this_an$datetime < first_depth_time))
          lidx <- min(which(this_an$datetime > last_depth_time))
          
          this_sub_anf <- this_an[c(fidx, fidx + 1), ] %>%
            dplyr::mutate(depth_min = depth - d_range,
                          depth_max = depth + d_range) %>%
            dplyr::select(deployid, datetime, depth_min, depth_max) %>%
            dplyr::summarise(depth_min = min(depth_min),
                             depth_max = max(depth_max))
          f_series_vec <- data.frame(deployid = animal, cee_id = cee_id, 
                                     datetime = first_depth_time,
                                     depth_min = this_sub_anf$depth_min,
                                     depth_max = this_sub_anf$depth_max)
          
          this_sub_anl <- this_an[c(lidx - 1, lidx), ] %>%
            dplyr::mutate(depth_min = depth - d_range,
                          depth_max = depth + d_range) %>%
            dplyr::select(deployid, datetime, depth_min, depth_max) %>%
            dplyr::summarise(depth_min = min(depth_min),
                             depth_max = max(depth_max))
          l_series_vec <- data.frame(deployid = animal, cee_id = cee_id, 
                                     datetime = last_depth_time,
                                     depth_min = this_sub_anl$depth_min,
                                     depth_max = this_sub_anl$depth_max)
          
          # Checking for duplication as some of the first or last times are the same
          dup_time_first <- f_series_vec$datetime == lubridate::round_date(f_series_vec$datetime, '5 minutes')
          dup_time_last <- l_series_vec$datetime == lubridate::round_date(l_series_vec$datetime, '5 minutes')
          if (dup_time_first & dup_time_last){
            
            this_an_series <- dplyr::bind_rows(this_an_series) %>%
              as.data.frame()
            
          } else if(!dup_time_first & dup_time_last){
            
            this_an_series <- dplyr::bind_rows(f_series_vec, this_an_series) %>%
              as.data.frame()
            
          } else {
            
            this_an_series <- dplyr::bind_rows(f_series_vec, this_an_series, l_series_vec) %>%
              as.data.frame()            
          }
          
        } # end if/else for non 2019 animals
        
        my_depths[[i]] <- this_an_series
        
      } # end loop over deploy ids
     
      an_series_out <- do.call("rbind", my_depths) %>% 
        arrange(deployid, cee_id, datetime)
      
      saveRDS(an_series_out, file = here::here('data/proc_data', 
                                           paste0('all_animals_depth_values.rds')))
      
      an_series_out
      
      
    })
)

