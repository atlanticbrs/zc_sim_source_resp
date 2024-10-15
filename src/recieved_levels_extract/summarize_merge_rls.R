register_target(
  
  # Here we do the merge of the Received Levels
  # calculated with the crawl-based rejection sampler
  # and also without. If there's ever a duplicate, I'll
  # keep the rejection sampler version
  
  tar_target(
    
    name = summarize_merge_rls,
    command = {
      
      # Let's bind the files
      # Rejection Sampler
      file_path <- list.files(here::here('data/proc_data'), 
                              pattern="*_rl_summary_rj.rds", full.names=TRUE)
      
      ## counting how many files in your directory
      file_number <- 1 
      
      for (file in file_path) {
        name <- paste0("df_", file_number) ## getting a new name ready
        x <- readRDS(file) ## reading in the file
        assign(name, x) ## assigning the data the name we created two lines up
        file_number <- file_number + 1 
      }
      
      full_data <- mget(ls(pattern = "^df_")) %>%
        reduce(bind_rows)
      
      rm(list = ls(pattern = "^df_"))
      
      df_sum <- full_data %>% 
        dplyr::filter(Quantile == 0.5) %>% 
        dplyr::group_by(deployid, cee_id) %>% 
        dplyr::slice(which.max(Value)) %>% 
        dplyr::mutate(max_val = TRUE)
      
      df_all <- left_join(full_data, df_sum,
                          by = c("deployid", "datetime")) %>% 
        dplyr::select(deployid, datetime, quantile = Quantile.x,
                      value = Value.x, cee_id = cee_id.x, max_val)
      
      df_plot <- df_all %>% 
        dplyr::filter(max_val == TRUE)
      
      df_pivot <- df_plot %>% 
        pivot_wider(
          names_from = quantile,
          values_from = value
        ) %>% 
        rename(median_rl = `0.5`, lower = `0.025`, upper = `0.975`)
      
      df_pivot$rl_bin <- cut(df_pivot$median_rl, breaks = c(0, 100, 120, 200), right = TRUE,
                             labels = c("<100", "100-120", ">120"))
      
      pivot_rj_df <- df_pivot %>% 
        dplyr::mutate(match_id = paste0(deployid, '_cee-id_', cee_id))
      rm(list = ls(pattern = "^df_"))
      
      # Regular
      file_path <- list.files(here::here('data/proc_data'), 
                              pattern="*_rl_summary.rds", full.names=TRUE)
      
      ## counting how many files in your directory
      file_number <- 1 
      
      for (file in file_path) {
        name <- paste0("df_", file_number) ## getting a new name ready
        x <- readRDS(file) ## reading in the file
        assign(name, x) ## assigning the data the name we created two lines up
        file_number <- file_number + 1 
      }
      
      full_data <- mget(ls(pattern = "^df_")) %>%
        reduce(bind_rows) %>% 
        drop_na(Value)
      
      rm(list = ls(pattern = "^df_"))
      
      df_sum <- full_data %>% 
        dplyr::filter(Quantile == 0.5) %>% 
        dplyr::group_by(deployid, cee_id) %>% 
        dplyr::slice(which.max(Value)) %>% 
        dplyr::mutate(max_val = TRUE)
      
      df_all <- left_join(full_data, df_sum,
                          by = c("deployid", "datetime")) %>% 
        dplyr::select(deployid, datetime, quantile = Quantile.x,
                      value = Value.x, cee_id = cee_id.x, max_val)
      
      df_plot <- df_all %>% 
        dplyr::filter(max_val == TRUE)
      
      df_pivot <- df_plot %>% 
        pivot_wider(
          names_from = quantile,
          values_from = value
        ) %>% 
        rename(median_rl = `0.5`, lower = `0.025`, upper = `0.975`)
      
      df_pivot$rl_bin <- cut(df_pivot$median_rl, breaks = c(0, 100, 120, 200), right = TRUE,
                             labels = c("<100", "100-120", ">120"))
      
      pivot_df <- df_pivot %>% 
        mutate(match_id = paste0(deployid, '_cee-id_', cee_id))
      
      # Now we merge
      sumidx <- !(pivot_df$match_id %in% pivot_rj_df$match_id)
      pivot_df <- pivot_df[sumidx, ]
      
      all_rls <- bind_rows(pivot_rj_df, pivot_df) %>% 
        dplyr::select(-match_id)
      
      f = file.path(here::here('results/output_data'))
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      write_csv(all_rls, file = file.path(f, paste0(Sys.Date(),'_rls_combined.csv')))
      
      all_rls
      
    })
  
  
)