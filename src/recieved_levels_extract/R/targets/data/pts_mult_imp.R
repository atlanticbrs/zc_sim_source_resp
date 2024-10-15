register_target(
  # 
  tar_target(
    name = pts_mult_imp,
    command = {

      
      # Below is the coding logic, but you can't call it like this
      # you need to use the manifest and loop one animal/cee combo
      # at a time. The code below will work when cee_times is constant
      # but yours needs to change with each animal/cee combo
      # list(samples = crawlUtils::cu_crw_sample(modelfits$fit,
      #                                          predTime = cee_times,
      #                                          size = 100),
      #      deployid = unique(modelfits$deployid))
      #
      # You need to work with these targets:
      # 1. cee_dat - data frame of metadata for each cee
      #
      # 2. cee_times - list of length 13
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
      
      # get the names of the model fit by peeling off the names of the sublist
      targetnames <- names(modelfits)
      animal_vec <- numeric()
      for (targetname in targetnames) {
        animal_vec <- c(animal_vec, modelfits[[targetname]][["deployid"]])
      }
      
      # pre-allocate the list
      mylist <- vector(mode = 'list', length = nrow(cee_dat))
      
      # convert the tibble to a data frame
      cee_df <- data.frame(cee_dat)
      
      for(i in 1:nrow(cee_df)){
        
        #pull out the animal and the cee
        animal <- cee_df[i, 'deployid'] 
        cee_id <- cee_df[i, 'cee_id']
        ser_status <- cee_df[i, 'ser']
        ser_gaps <- cee_df[i, 'ser_gaps']
        
        # find the times for that cee
        my_cee <- which(cee_id == names(cee_times))
        
        # get the times
        my_cee_times <- cee_times[[my_cee]]$times
        
        # Model fit extraction
        my_model_fit_id <- which(animal == animal_vec)
        my_model_fit <- modelfits[[my_model_fit_id]][["fit"]]
      
        # make the predictions - 100 by the length of my_cee_times as well as observed points
        samples = crawlUtils::cu_crw_sample(my_model_fit, 
                                            predTime = my_cee_times,
                                            size = 100)
        
        # Make samples into a smaller data frame
        mydf <- dplyr::bind_rows(samples) %>% 
          dplyr::filter(locType == 'p') %>% 
          dplyr::mutate(deployid = animal,
                 cee_id = cee_id,
                 ser_status = ser_status,
                 ser_gaps = ser_gaps)
        
        # return the values
        mylist[[i]] <- list(deployid = animal,
                            cee_id = cee_id,
                            series_stats = ser_status,
                            series_gaps = ser_gaps,
                            pred_times = my_cee_times,
                            imputed_tracks = samples,
                            predicted_pts = mydf)
          
      }
      
      saveRDS(mylist, file = here::here('data/proc_data', 'imputed_tracks.rds'))
      
      mylist
      
    }
  )
)