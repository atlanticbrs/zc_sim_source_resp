register_target(
  # 
  tar_target(
    name = run_hanks_short,
    command = {
      
      # Goal is to run the Hanks model
      # depends on the hanks_data_prep target
      
      # pre-allocate the list
      mylist <- vector(mode = 'list', length = length(hanks_data_prep_short))
      pweight <- 32
      
      for(i in 1:length(mylist)){
        
        my_deployid <- hanks_data_prep_short[[i]][["deployid"]]
        my_cee_id <- hanks_data_prep_short[[i]][["cee_id"]]
        my_glm_data <- hanks_data_prep_short[[i]][["hanks_glm_data"]]
        
        #assemble a data frame
        # Get correct covariates out
        d2ship_idx <- which(stringr::str_detect(names(my_glm_data), 
                                                glob2rx("d2ship*grad")))
        shlf_idx <- which(names(my_glm_data) == "shelf_brk_loc")
        gam_dat_48 <- my_glm_data %>%
          dplyr::select(
            c(z,
              tau,
              t,
              datetime,
              time_scl,
              d2ship_grad = names(my_glm_data)[d2ship_idx],
              shelf_break_loc = names(my_glm_data)[shlf_idx],
              crw,
              exposure,
              exposure_ordered,
              t_exposure
            )
          ) %>%
          dplyr::mutate(d2ship_grad_exp = ifelse(exposure == 'YES', d2ship_grad, 0))%>% 
          data.frame()
        
        # Initialize an empty list for the current animal
        current_results <- list(deployid = my_deployid, cee_id = my_cee_id)
        
        tick <- Sys.time()
        tryCatch({
          # First gam fitting
          # Varying coefficient for dist 2 ship | t
          gam_result_1 <- tryCatch({
            fit_1 = gam(
              z ~ s(t, by = d2ship_grad, k = 40, bs = "ad") + 
                crw,
              weights = rep(1 / pweight, nrow(gam_dat_48)),
              family = "poisson",
              # method = 'REML',
              offset = log(tau),
              data = gam_dat_48
            )
          }, error = function(e) {
            cat("Error in GAM 1 for subset: ", my_deployid, "\n")
            cat("Error message: ", e$message, "\n")
            NULL
          })
          
          if (!is.null(gam_result_1)) {
            saveRDS(fit_1, file = here::here(
              'results/gams/cee_48',
              paste0(my_deployid, "_",
                     my_cee_id,
                     '_gam_1.rds')
            ))
          }
          
         # Add the GAM results to the current animal's results
          current_results$fit_1 <- gam_result_1
         
        }, finally = {
          tock <- Sys.time() - tick
          duration_minutes <- as.numeric(tock, units = "mins")
          current_results$time_to_fit <- duration_minutes
          
          # return the values
          mylist[[i]] <- current_results
          
        }) 
      }
      
      return(mylist)
    }
    
  )
  
)