register_target(
  # UPDATE this HELP
  tar_target(
    name = locs_data,
    command = {
      
      my_cols <- cols(
        DeployID = col_character(),
        Ptt = col_double(),
        Instr = col_character(),
        Date = col_double(),
        Type = col_character(),
        Quality = col_character(),
        Latitude = col_double(),
        Longitude = col_double(),
        Error.radius = col_double(),
        Error.Semi.major.axis = col_double(),
        Error.Semi.minor.axis = col_double(),
        Error.Ellipse.orientation = col_double(),
        Offset = col_logical(),
        Offset.orientation = col_logical(),
        GPE.MSD = col_logical(),
        GPE.U = col_logical(),
        Count = col_logical(),
        Comment = col_character(),
        originaldate = col_character()
      )
      
      f = file.path(here::here('data'))
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      tbl_locs <- list.files(here::here("../../01_shared_data_products/filter_sattag/locations"),
                             recursive = TRUE, 
                             full.names = TRUE, 
                             pattern = "*-Locations.csv") %>% 
        purrr::map_df( ~ readr::read_csv(.x, col_types = my_cols)) %>% 
        janitor::clean_names() %>%
        dplyr::rename(deployid = deploy_id) %>% 
        dplyr::mutate(datetime = as.POSIXct(date, tz = "UTC", origin = "1970-01-01")) %>%
        readr::write_csv(file = file.path(f, "ss_locs.csv"))
      
      
      
      tbl_locs
      
    }
  )
)