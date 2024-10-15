register_target(
  # 
  tar_target(
    name = hanks_data_prep_short,
    command = {
      
      # Goal is to prep the data for Hanks by taking imputations
      # from the pts_mult_imp_hanks_short target and doing the processing
      # to convert them to the ctmc

      # pre-allocate the list
      mylist <- vector(mode = 'list', length = length(pts_mult_imp_hanks_short))
      
      for(i in 1:length(mylist)){
        
        #pull out the animal and the cee
        animal <- pts_mult_imp_hanks_short[[i]][["deployid"]]
        my_cee_id <- pts_mult_imp_hanks_short[[i]][["cee_id"]]
        my_raster_name <- paste0('d2ship_cee_', my_cee_id)
        pred_pts <- pts_mult_imp_hanks_short[[i]][["predicted_pts"]]
        
        print(paste(i, animal, my_cee_id, sep = ": "))
        
        # Set up a set of time ranges to define exposure periods, and 
        # error trap for multiple cees
        cee_start_time <- cee_dat %>% 
          dplyr::filter(deployid == animal,
                        cee_id == my_cee_id) %>% 
          dplyr::mutate(exp_window = start_datetime + lubridate::hours(24))
        
        # Here's where we'd cull the data somehow based on what we set up as our exposure window
        # First we set up the temporal extent of it
        cee_48 <- data.frame(cee_window_start = cee_start_time$start_datetime - lubridate::hours(24), 
                             cee_window_end = cee_start_time$start_datetime + lubridate::hours(24))
        
        if(cee_start_time$start_datetime - range(pred_pts$datetime)[1] < lubridate::hours(24)) next()
        
        # Need to drop the geometry now because of a conflict with the vectors library and sf objects
        geometry <- sf::st_geometry(pred_pts)
        pred_pts_df <- sf::st_set_geometry(pred_pts, NULL)
        idx <- which(pred_pts_df$datetime >= cee_48$cee_window_start & pred_pts_df$datetime < cee_48$cee_window_end)
        
        pred_pts_df <- pred_pts_df %>%
          dplyr::filter(datetime >= cee_48$cee_window_start & datetime < cee_48$cee_window_end)
        
        pred_pts <- sf::st_set_geometry(pred_pts_df, geometry[idx])
        
        cee_exp_window <- cee_start_time[, c('start_datetime', 'exp_window', 'cee_id', 'cee_type'), drop = TRUE] %>% 
          dplyr::mutate(focal_cee = if_else(cee_id == my_cee_id, "YES", "NO"))

        # pare down to just one imputation
        glm_data <- NULL
        for(j in unique(pred_pts$rep)){
          
          my_pts_sub <- pred_pts %>% 
            filter(rep == j)
          
          # run the Hanks Prep
          crop_lim <- raster::extent(c(range(sf::st_coordinates(my_pts_sub)[, 'X']) + c(-5000, 5000), 
                                       range(sf::st_coordinates(my_pts_sub)[, 'Y']) + c(-5000, 5000)))
          grad.stack_crop <- raster::stack(raster::crop(prep_rasters[[which(names(prep_rasters) == my_raster_name)]], crop_lim),
                                           raster::crop(prep_rasters[['bathy']], crop_lim)) # raster stack of gradient covariates
          names(grad.stack_crop) <- paste0(names(grad.stack_crop), "_grad")
          loc.stack_crop <- raster::stack(raster::crop(prep_rasters[['shelf_brk']], crop_lim)) # raster stack of motility covariates                      
          names(loc.stack_crop) <- paste0(names(loc.stack_crop), "_loc")
          
          path <- list(xy = sf::st_coordinates(my_pts_sub), t = as.vector(my_pts_sub$datetime))
          ctmc <- ctmcmove::path2ctmc(xy = path$xy, t = path$t, rast = grad.stack_crop) # extract discrete path & cell residence times
          glm_data <- rbind(glm_data, ctmcmove::ctmc2glm(ctmc, loc.stack_crop, grad.stack_crop)) 
          
        }
        
        # Establish the baseline vs exposure factor covariate
        glm_data$datetime <- as.POSIXct(glm_data$t, origin = '1970-01-01', tz = 'UTC')
        
        # Establish the Before/After Exposure Factor, an ordered version, and a time since variable
        # Also add a rescaled time so I can plot all the time varying coefficients on top of each other
        glm_data <- glm_data %>% 
          dplyr::mutate(exposure = forcats::as_factor(if_else(datetime > cee_start_time$start_datetime, 'YES', 'NO')),
                        t_exposure = 0,
                        time_scl = t - min(t))
        glm_data$exposure_ordered <- factor(glm_data$exposure, ordered = TRUE)
        
        # Calculate the elapsed time
        glm_data <- glm_data %>% 
          dplyr::mutate(reset_accumulator = cumsum(exposure == 'YES' & lag(exposure, default = 'NO') == 'NO')) %>%
          dplyr::mutate(t_exposure = ifelse(exposure == 'YES', 
                                            cumsum(as.numeric(difftime(datetime, lag(datetime, 
                                                                                    default = first(datetime)), 
                                                                       units = "hours"))), 0))  %>% 
          dplyr::ungroup() %>%
          dplyr::select(-reset_accumulator) 
        
        # return the values
        mylist[[i]] <- list(deployid = animal,
                            cee_id = my_cee_id,
                            hanks_glm_data = glm_data)
        
      }
      proc_list <- vctrs::list_drop_empty(mylist)
      return(proc_list)
      
    }
  )
)