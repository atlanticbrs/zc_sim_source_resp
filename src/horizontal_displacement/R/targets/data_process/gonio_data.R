register_target(
  # set value = TRUE to build targets without updating individual targets 
  # that are expensive to run, such as MCMC samplers.  such targets must be 
  # designed to return non-breaking output before expensive steps are to be run
  #  if skip_expensive_targets$value == TRUE.
  tar_target(
    name = gonio_data,
    command = {
      
      # Listing out the files since there's a change in the DeployID to deployid column that messes up
      # map_df
      gonio_2017_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2017.csv")
      gonio_2018_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2018.csv")
      gonio_2019_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2019.csv")
      gonio_2020_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2020.csv")     
      gonio_2021_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2021.csv")
      gonio_2022_file <-  here::here("../../01_shared_data_products/gonio/gonio_gpx_merge_2022.csv")
      
      # read the files in, in turn
      gonio_dat_2017 <- read_csv(gonio_2017_file, show_col_types = FALSE) %>% 
        janitor::clean_names() %>% 
        dplyr::rename(deployid = deploy_id)
      
      gonio_dat_2018 <- read_csv(gonio_2018_file, show_col_types = FALSE) %>% 
        janitor::clean_names() %>% 
        dplyr::rename(deployid = deploy_id)
      
      gonio_dat_2019 <- read_csv(gonio_2019_file, show_col_types = FALSE) %>% 
        janitor::clean_names() %>% 
        dplyr::rename(deployid = deploy_id)
      
      gonio_dat_2020 <- read_csv(gonio_2020_file, show_col_types = FALSE) %>% 
        janitor::clean_names() %>% 
        dplyr::rename(deployid = deploy_id)
      
      gonio_dat_2021 <- read_csv(gonio_2021_file, show_col_types = FALSE) %>% 
        janitor::clean_names() 
      
      gonio_dat_2022 <- read_csv(gonio_2022_file, show_col_types = FALSE) %>% 
        janitor::clean_names() 
      
      # Bind together
      gonio_dat <- dplyr::bind_rows(gonio_dat_2017, gonio_dat_2018, 
                                    gonio_dat_2019, gonio_dat_2020, 
                                    gonio_dat_2021, gonio_dat_2022) %>% 
        dplyr::mutate(deployid = stringr::str_remove(deployid, "_DUML")) %>% 
        dplyr::filter(deployid %in% exposed_ans) %>% 
        dplyr::mutate(strength_db = if_else(strength_db == 0, -30, strength_db)) %>% 
        dplyr::filter(strength_db != 1)
      
      readr::write_csv(gonio_dat, file = here::here('data',
                                           'gonio_locs.csv'))
      
      gonio_dat

      
    }
  )
)