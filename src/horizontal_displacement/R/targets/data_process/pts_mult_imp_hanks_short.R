register_target(
  # 
  tar_target(
    name = pts_mult_imp_hanks_short,
    command = {

      # Goal is step 1 of the data prep for Hanks
      # Specifically, this will create 32 imputations 
      # at a 15 minute temporal resolution
      # 2024-02-21: this is for the 24/24 test for Hanks 
      
      # get the names of the model fit by peeling off the names of the sublist
      targetnames <- names(modelfits)
      animal_vec <- numeric()
      for (targetname in targetnames) {
        animal_vec <- c(animal_vec, modelfits[[targetname]][["deployid"]])
      }
      
      # pre-allocate the list
      mylist <- vector(mode = 'list', length = length(modelfits))
      
      # convert the tibble to a data frame
      cee_df <- data.frame(cee_dat)
      
      for(i in 1:nrow(cee_df)){
        
        #pull out the animal and the cee
        animal <- cee_df[i, 'deployid'] 
        cee_id <- cee_df[i, 'cee_id']
        ser_status <- cee_df[i, 'ser']
        ser_gaps <- cee_df[i, 'ser_gaps']
        
        # Model fit extraction
        my_model_fit_id <- which(animal == animal_vec)
        my_model_fit <- modelfits[[my_model_fit_id]][["fit"]]
      
        # make the predictions - 32 by the length of my_cee_times as well as observed points
        samples = crawlUtils::cu_crw_sample(my_model_fit, 
                                            predTime = "15 min",
                                            size = 32)
        
        # Make samples into a smaller data frame
        mydf <- dplyr::bind_rows(samples) %>% 
          dplyr::mutate(deployid = animal,
                 cee_id = cee_id,
                 ser_status = ser_status,
                 ser_gaps = ser_gaps)
        
        # return the values
        mylist[[i]] <- list(deployid = animal,
                            cee_id = cee_id,
                            predicted_pts = mydf)
          
      }
      
      return(mylist)
      
    }
  )
)