register_target(
  
  # Here we do the summary of the RL data _after_
  # we've done the sampling for 
  
  tar_target(
    
    name = summarize_ambient,
    command = {
      
      
      # Container to hold output
      rl_summary <- vector(mode = 'list', length = length(sample_ambient))
      # loop over one animal/cee combination at a time
      for(i in 1:length(sample_ambient)){
        
        my_cee_id <- sample_ambient[[i]][["cee_id"]]
        my_deployid <- sample_ambient[[i]][["deployid"]]
        my_cee_type <- cee_dat %>% 
          dplyr::filter(cee_id == my_cee_id,
                        deployid == my_deployid) %>% 
          dplyr::select(cee_type)
        # We can't summarize controls - no RLs
        if(as.logical(my_cee_type == 'control')) next()
        animal <- sample_ambient[[i]][["pts_tobl"]]
        
        print(paste(i, my_deployid, my_cee_id, sep = " - "))
        
        # Do the Re-sampling on RL
        rl_boot_samples <- vector(mode = 'list', length = nrow(animal))
        rl_samples <- vector(mode = 'list', length = nrow(animal))
        for(j in 1:nrow(animal)){
          
          # Grab the vector of RLs through the H20 col at a given x,y,t
          rls <- unlist(animal$all_rl_vals[j])
          if(length(rls) == 0) next()
          # Subtract the 1/3 Octave band 
          # if negative, than source was not detectable (NA)
          # if positive, retain
          rl_amb <- rls - animal$tobl[j]
          rls <- if_else(rl_amb <= 0, NA, rls)
          
          # if all values are NA, go to the next row
          if(all(is.na(rl_amb))) {
            rl_smpl <- NA
          } else {
            # sample with replacement from rl vector & write out to a data frame
            rl_smpl <- sample(rls, 500, replace = TRUE)  
          }
          
          # sample with replacement from rl vector
          rl_samples[[j]] <- data.frame(deployid = animal$deployid[j],
                                        datetime = animal$datetime[j],
                                        cee_id = animal$cee_id[j],
                                        track_id = animal$rep[j],
                                        raw = rls)
          
          rl_boot_samples[[j]] <- data.frame(deployid = animal$deployid[j],
                                        datetime = animal$datetime[j],
                                        cee_id = animal$cee_id[j],
                                        track_id = animal$rep[j],
                                        raw = rl_smpl)
          
        }
        rl_df <- do.call('rbind', rl_samples)
        rl_boot_df <- do.call('rbind', rl_boot_samples)
        
        f = file.path(here::here('data/proc_data/rl_samples'))
        dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
        
        write_csv(rl_df, file = file.path(f, paste0(my_deployid, "_", my_cee_id, "_rl-samples.csv")))
                  
        # End Re-sampling on RL
          
        # Write the RL Quantile Summary to Disk
        rl_sum <- rl_boot_df %>% 
          group_by(deployid, datetime) %>% 
          do({x <- .$raw
          map_dfr(.x = c(.025, .5, .975),
                  .f = ~ tibble(Quantile = .x,
                                    Value = quantile(x, probs = .x, na.rm = TRUE)))
          }) %>% 
          mutate(cee_id = my_cee_id)
        saveRDS(rl_sum, file = here::here('data/proc_data', 
                                         paste0(unique(animal$deployid), "_", 
                                                my_cee_id, '_rl_summary.rds')))
        
        # Quantile by track ID for cSEL cals
        rl_sum_track <- rl_df %>% 
          group_by(deployid, datetime, track_id) %>% 
          do({x <- .$raw
          map_dfr(.x = c(.025, .5, .975),
                  .f = ~ tibble(Quantile = .x,
                                Value = quantile(x, probs = .x, na.rm = TRUE)))
          }) %>% 
          mutate(cee_id = my_cee_id)
        
        # Make and save an exploratory plot
        p <- ggplot2::ggplot(rl_boot_df, aes(x = raw, y = factor(datetime), 
                                             fill = after_stat(x)))+
          ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
          theme_bw()+
          labs(title = "Estimated RLs at 5' Intervals within the CEE",
               subtitle = paste0(unique(animal$deployid), ", CEE ID: ", my_cee_id),
               y = '', x = expression(paste("SPL (dB re:1 ", mu, "Pa)")),
               fill = 'SPL')+
        scale_fill_viridis_c()
        
        f = file.path(here::here('results/plots/rl_summary'))
        dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
        
        ggsave(p, 
               filename = file.path(f, paste0(unique(animal$deployid), "_", my_cee_id, "_rl-summary.png")),
               device = 'png',
               dpi = 'retina',
               width = 10, height = 6.1, units = 'in') 
        # End Plotting
        
        
        # Store items in List 
        rl_summary[[i]] <- list(deployid = unique(animal$deployid),
                            cee_id = my_cee_id,
                            rl_raw = animal,
                            rl_summary = rl_df,
                            rl_quantiles = rl_sum,
                            rl_quantiles_bytrack = rl_sum_track)
        
      } # end i loop over the list of each animal/cee combination
      
      rl_summary
      
      })
             

)