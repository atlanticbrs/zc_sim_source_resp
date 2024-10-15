register_target(
  # 
  tar_target(
    name = pts_mult_imp_rejsamp,
    command = {
      
      # Load in the bathy layer
      new_bathy_prj  <- readRDS(file = bathy_file_location)

      # You need to work with these targets:
      # 1. cee_dat - data frame of metadata for each cee, but we can only use series data for these
      #
      # 2. cee_times - list of length animals
      #     a) cee_times[[1]][[1]]
      #     [1] "17_01"
      #     b) > cee_times[[1]][["times"]]
      #     [1] "2017-08-22 18:41:00 UTC" "2017-08-22 18:40:00 UTC" "2017-08-22 18:45:00 UTC" "2017-08-22 18:50:00 UTC"
      #     [5] "2017-08-22 18:51:00 UTC"
      #
      # 3. modelfits - list of length 50 (animals); each sublist is length 2:
      #     a) modelfits[["modelfits_6ed74ead"]][["fit"]]
      #         <...>
      #     b) modelfits[["modelfits_6ed74ead"]][["deployid"]]
      #       [1] "ZcTag060"
      # 4. series_data - contains tbl_series which has the depth information
      #     from all animals. Will use this for 
      #### Processing
      #
      # Note that unlike pts_mult_imp, here we will employ the rejection sampler, so we need to get the 
      #     depth data matched up to the times
      
      # get the names of the model fit by peeling off the names of the sublist
      targetnames <- names(modelfits)
      animal_vec <- numeric()
      for (targetname in targetnames) {
        animal_vec <- c(animal_vec, modelfits[[targetname]][["deployid"]])
      }
      
      # Only use animals that have series data AND
      # only use animals whose series data is active at time of CEE AND
      # which don't have gaps at the cee
      # TODO: should cee_rj be its own target?
      cee_rj <- cee_dat %>% 
        filter(ser, 
               cee_type == 'simulated mfas',
               cee_st < ser_en,
               is.na(ser_gaps))
      
      # pre-allocate the list to contain output for each animal
      mylist <- vector(mode = 'list', length = nrow(cee_rj))
      
      # convert the tibble to a data frame
      cee_df <- data.frame(cee_rj)
      
      for(i in 1:nrow(cee_df)){
        
        #pull out the animal and the cee
        animal <- cee_df[i, 'deployid'] 
        my_cee_id <- cee_df[i, 'cee_id']
        
        print(paste(animal, my_cee_id, sep = ' - '))
        
        # find the times for that cee
        my_cee <- which(my_cee_id == names(cee_times))
        
        # get the times
        my_cee_times <- cee_times[[my_cee]]$times
        
        # Model fit extraction
        my_model_fit_id <- which(animal == animal_vec)
        my_model_fit <- modelfits[[my_model_fit_id]][["fit"]]
        
        # Get the Depth Data
        an_series <- an_series_df %>% 
          filter(deployid == animal,
                 cee_id == my_cee_id)
        
        #######################################################
        # make the predictions with the rejection sampler     #
        # This list (outlist) will contain 100 predicted points at each
        # of the time points within a cee.
        outlist <- vector(mode = 'list', length = nrow(an_series))
        for(j in 1:nrow(an_series)){
          
          # List to contain 100 locations for each animal
          n_imps <- 100
          imputes <- vector(mode = 'list', length = n_imps) 
          
          n_impute <- 1
          while(n_impute <= n_imps){ 
            
            # Simulate new points from the fitted crawl object
            # New way of imputing
            crw_sim = crawlUtils::cu_crw_sample(my_model_fit, 
                                                predTime = an_series[j, 'datetime'],
                                                size = 1) 
            
            new_pt <- crw_sim[[1]] %>%
              dplyr::filter(locType == "p") %>% 
              dplyr::mutate(imputation_num = paste0('track_', n_impute),
                            bathymetry = NA) 
            
            # Sample bathy raster, and test for depth mismatch
            new_pt$bathymetry <- raster::extract(new_bathy_prj, new_pt) 
            df_join <- left_join(an_series[j, ], new_pt, by = 'datetime')
            
            df <- df_join %>% 
              dplyr::mutate(keep = depth_max < abs(bathymetry))
            
            
            # append good rows, i.e. those that are all in good depth ranges
            if(all(df$keep == TRUE, na.rm = TRUE)) {
              imputes[[n_impute]] <- df
              # imputes_df <- bind_rows(imputes_df, df)
              n_impute <- n_impute + 1
              # print(n_impute)
            } else {
              next
            }
            
          } # end while loop over the imputations
          
          # put the data frame of the 100 "kept" points into a list
          # recall this list is the length of the predicted
          # time vector
          imputes_df <- do.call("rbind", imputes) 
          outlist[[j]] <- imputes_df 
          
        } # end j loop over animal times
        df <- do.call("rbind", outlist) 
        
        # End rejection sampler                               #
        #######################################################
        
        # return the values
        mylist[[i]] <- list(deployid = animal,
                            cee_id = my_cee_id,
                            an_depth_data = an_series,
                            pred_times = my_cee_times,
                            imputed_points = df)
          
      }
      
      saveRDS(mylist, file = here::here('data/proc_data', 'imputed_tracks_rj.rds'))
      return(mylist)
      
    }
  )
)