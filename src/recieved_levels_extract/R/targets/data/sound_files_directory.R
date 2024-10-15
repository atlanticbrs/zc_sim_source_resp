register_target(
  
  # Locations for:
  # cee file that has the water column depths
  # for all the animals in 2019
  
  tar_target(
    
    name = sound_files_directory,
    command = {
      
      cee_1701_rl <- dir(here::here("../../02_large_data/cee_1701"), 
                      pattern = "^CEE.*\\.nc$",
                      recursive = TRUE,
                      full.names = TRUE)
      
      cee_1802_rl <- dir(here::here("../../02_large_data/cee_1802"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1803_rl <- dir(here::here("../../02_large_data/cee_1803"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1808_rl <- dir(here::here("../../02_large_data/cee_1808"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1901_rl <- dir(here::here("../../02_large_data/cee_1901"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1902_rl <- dir(here::here("../../02_large_data/cee_1902"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1903_rl <- dir(here::here("../../02_large_data/cee_1903"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_1904_rl <- dir(here::here("../../02_large_data/cee_1904"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      cee_2003_rl <- dir(here::here("../../02_large_data/cee_2003"),
                         pattern = "^CEE.*\\.nc$",
                         recursive = TRUE,
                         full.names = TRUE)
      
      rl_file_list <- list(list(cee_id = '17_01', 
                                rl_cee_times = cee_times$`17_01`$times,
                                cee_1701_files = cee_1701_rl),
                           list(cee_id = '18_02', 
                                rl_cee_times = cee_times$`18_02`$times,
                                cee_1802_files = cee_1802_rl),
                           list(cee_id = '18_03', 
                                rl_cee_times = cee_times$`18_03`$times,
                                cee_1803_files = cee_1803_rl),
                           list(cee_id = '18_08', 
                                rl_cee_times = cee_times$`18_08`$times,
                                cee_1808_files = cee_1808_rl),
                           list(cee_id = '19_01', 
                                rl_cee_times = cee_times$`19_01`$times,
                                cee_1901_files = cee_1901_rl),
                           list(cee_id = '19_02', 
                                rl_cee_times = cee_times$`19_02`$times,
                                cee_1902_files = cee_1902_rl),
                           list(cee_id = '19_03', 
                                rl_cee_times = cee_times$`19_03`$times,
                                cee_1903_files = cee_1903_rl),
                           list(cee_id = '19_04', 
                                rl_cee_times = cee_times$`19_04`$times,
                                cee_1904_files = cee_1904_rl),
                           list(cee_id = '20_03', 
                                rl_cee_times = cee_times$`20_03`$times,
                                cee_2003_files = cee_2003_rl))
      
      rl_file_list
      
      })
             

)