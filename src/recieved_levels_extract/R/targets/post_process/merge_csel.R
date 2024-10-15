register_target(
  
  # Here we do the merge of the cSEL files
  
  tar_target(
    
    name = merge_csel,
    command = {
      
      # Let's bind the files
      # Rejection Sampler

      calc_cSEL_rj_df <- calc_cSEL_rj %>% 
        mutate(match_id = paste0(deployid, cee_id))
      
      calc_cSEL_df <- calc_cSEL %>% 
        mutate(match_id = paste0(deployid, cee_id))
      
      # Now we merge - logic is I want the ones from the _non_ rejection-sampler
      # data frame that are not also in the rj data frame. This means that if an
      # animal is in both, I keep the one in the rj data frame
      sumidx <- !(calc_cSEL_df$match_id %in% calc_cSEL_rj_df$match_id)
      calc_cSEL_df <- calc_cSEL_df[sumidx, ]
      
      all_csel <- bind_rows(calc_cSEL_rj_df, calc_cSEL_df) %>% 
        dplyr::select(-match_id) %>% 
        rename(lower = `0.025`, median = `0.5`, upper = `0.975`)
      
      write_csv(all_csel, here::here('results/output_data/cSEL',
                                     paste0(Sys.Date(),'_csel_combined.csv')))
      

      all_csel
    })
  
  
)