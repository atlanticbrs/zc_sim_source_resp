register_target(
  # THis is to read in the series data
  tar_target(
    name = series_data,
    command = {
      
      my_cols <- cols(
        DeployID = col_character(),
        Ptt = col_double(),
        DepthSensor = col_double(),
        Source = col_character(),
        Instr = col_character(),
        Day = col_character(),
        Time = col_time(format = ""),
        LocationQuality = col_logical(),
        Latitude = col_logical(),
        Longitude = col_logical(),
        Depth = col_double(),
        DRange = col_double(),
        Temperature = col_logical(),
        TRange = col_logical(),
        Activity = col_logical(),
        ARange = col_logical(),
        Date = col_double(),
        original = col_character()
      )
      
      tbl_series <- list.files(here::here("../../01_shared_data_products/filter_sattag/series"),
                             recursive = TRUE, 
                             full.names = TRUE, 
                             pattern = "*-Series.csv") %>% 
        purrr::map_df( ~ readr::read_csv(.x, col_types = my_cols)) %>% 
        janitor::clean_names() %>%
        dplyr::rename(deployid = deploy_id) %>% 
        dplyr::mutate(deployid = stringr::str_remove(deployid, "_DUML")) %>% 
        dplyr::mutate(datetime = as.POSIXct(date, tz = "UTC", origin = "1970-01-01")) %>%
        readr::write_csv(file = here::here('data/',
                                           'ss_series.csv'))
      
      tbl_series
      
    }
  )
)