register_target(
  # Ok, here we pass in the locs_data and 
  # simply find the unique animals therein
  # these are the exposed animals (here called
  # an R object named deployids)
  tar_target(
    name = exposed_ans,
    command = {
      
      data <- locs_data %>% 
        dplyr::filter(species == 'Zc') %>% 
        dplyr::filter(deployid %in% cee_dat$deployid)
      deployids <- unique(data$deployid)
        
      saveRDS(deployids, file = here::here('data', 'deployids.rds'))
      
      deployids
      
    }
  )
)