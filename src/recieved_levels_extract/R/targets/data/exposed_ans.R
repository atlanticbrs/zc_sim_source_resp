register_target(
  # Ok, here we pass in the locs_data and 
  # simply find the unique animals therein
  # these are the exposed animals (here called
  # an R object named deployids)
  tar_target(
    name = exposed_ans,
    command = {
      
      data <- locs_data # here we have to pass in the value 
      deployids <- stringr::str_remove(unique(data$deployid), "_DUML")
      saveRDS(deployids, file = here::here('data/', 'deployids.rds'))
      
      deployids
      
    }
  )
)