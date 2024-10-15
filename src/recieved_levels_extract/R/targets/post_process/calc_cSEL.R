register_target(
  
  # Here we assemble the cSEL vectors. For the non-rejection sampler
  
  tar_target(
    
    name = calc_cSEL,
    command = {
      
      compact_summarize_ambient <- compact(summarize_ambient)
      cSEL_summary <- vector(mode = "list", 
                             length = length(compact_summarize_ambient))
      for(i in 1:length(compact_summarize_ambient)){
        
        my_deployid <- compact_summarize_ambient[[i]][["deployid"]]
        my_cee_id <- compact_summarize_ambient[[i]][["cee_id"]]
        this_cee <- cee_times[[my_cee_id]][["times"]]
        sonar_seq <- seq.POSIXt(this_cee[1], 
                                this_cee[length(this_cee)], 
                                by = '25 s')
        
        # some of these are empty because they are controls
        rls <- compact_summarize_ambient[[i]][["rl_quantiles_bytrack"]] 
        
        if(nrow(rls) == 0) {
          print(paste(my_deployid, my_cee_id, sep = ': '))
          next()
        } else {
          rls <- rls %>%
            tidyr::drop_na(Value) %>% 
            dplyr::filter(Quantile == 0.5) %>% 
            dplyr::select(-Quantile) %>% 
            data.frame()
        }
        
        # some were really far and the above call to 
        # drop_na() zeros out the rls data frame
        if(nrow(rls) == 0) {
          print(paste(my_deployid, my_cee_id, sep = ': '))
          next()
        }
        
        csel_list<- vector(mode = 'list', length = nrow(rls))
        for(j in unique(rls$track_id)){
          rl_sub <- rls %>% 
            dplyr::filter(track_id == j)
          if(nrow(rl_sub) < 2) next()
          lin_int <- approx(rl_sub$datetime, 
                            rl_sub$Value, n = length(sonar_seq))
          
          csel_lin <- 10*log10(sum(10^((lin_int$y) / 10)))
          
          csel_list[[j]] <- data.frame(deployid = unique(rl_sub$deployid),
                                       track_id = unique(rl_sub$track_id),
                                       cee_id = unique(rl_sub$cee_id),
                                       cSEL = csel_lin)
          
        }
        
        csel_df <- do.call('rbind', csel_list)
        
        csel_sum <- csel_df %>% 
          do({x <- .$cSEL
          map_dfr(.x = c(.025, .5, .975),
                  .f = ~ tibble(Quantile = .x,
                                Value = quantile(x, probs = .x, na.rm = TRUE)))
          }) %>% 
          dplyr::mutate(deployid = my_deployid,
                        cee_id = my_cee_id)
        
        
        cSEL_summary[[i]] <- csel_sum
        
      }
      
      cSEL_summary <- do.call('rbind', cSEL_summary) %>% 
        pivot_wider(names_from = Quantile,
                    values_from = Value)
      
      f = file.path(here::here('results/output_data/cSEL'))
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      write_csv(cSEL_summary, file = file.path(f, "csel_summary.csv"))
      
      cSEL_summary
      
    })
  
)