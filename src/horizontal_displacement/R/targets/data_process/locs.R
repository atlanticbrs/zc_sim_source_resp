register_target(
  # Read in three targets and merge them together
  
  tar_target(
    name = locs,
    # command = wrangle_data(locs_data, gonio_data, deploy_file)
    command = {
      
      deploy_cols <- cols(
        ptt = col_factor(),
        hex = col_character(),
        serial = col_character(),
        tag_model = col_character(),
        deployid = col_character(),
        catalog_ID_DUML = col_character(),
        BRS_season = col_character(),
        deploydate_yyyymmdd = col_double(),
        deploytime_hhmm = col_double(),
        species = col_character(),
        location = col_character(),
        tagger = col_character(),
        permit = col_character(),
        tagtype = col_character(),
        longitude = col_double(),
        latitude = col_double(),
        tagside = col_character(),
        taglocation = col_character(),
        taglocationnotes = col_character(),
        programming_basic_descriptor = col_character(),
        notes = col_character(),
        recorder = col_character(),
        age = col_character(),
        sex = col_character(),
        agesex_conf = col_character(),
        headshot = col_character(),
        agesex_reason = col_character(),
        appearance_notes = col_character()
      )
      
      locs <- locs_data %>% 
        dplyr::filter(quality != 'Z') %>% 
        dplyr::mutate(ptt = forcats::as_factor(ptt),
                      year = lubridate::year(datetime),
                      month = lubridate::month(datetime),
                      day = lubridate::day(datetime),
                      quality = factor(quality,
                                       levels = c("3","2","1","0","A","B"))) %>%
        dplyr::mutate(
          error_semi_major_axis = ifelse(type == 'FastGPS', 50,
                                         error_semi_major_axis),
          error_semi_minor_axis = ifelse(type == 'FastGPS', 50,
                                         error_semi_minor_axis),
          error_ellipse_orientation = ifelse(type == 'FastGPS', 0,
                                             error_ellipse_orientation)
        ) %>% 
        dplyr::filter(data.table::between(latitude, 30, 42) & data.table::between(longitude, -77.5, -70))
      
      # Here we add the deployment location data
      deploy_dat <- readr::read_csv(deploy_file_location,
                                       col_types = deploy_cols,
                                    show_col_types = FALSE) %>% 
        dplyr::mutate(deployid = stringr::str_remove(deployid, "_DUML")) %>% 
        filter(deployid %in% unique(locs$deployid)) %>% 
        dplyr::mutate(longitude = ifelse(longitude > 0, -1 * longitude, longitude)) %>% 
        dplyr::mutate(datetime = lubridate::ymd_hm(paste0(deploydate_yyyymmdd, deploytime_hhmm))) %>% 
        dplyr::mutate(quality = 'User',
               type = 'User',
               error_semi_major_axis = 50,
               error_semi_minor_axis = 50,
               error_ellipse_orientation = 0) %>% 
        dplyr::select(c("deployid", "datetime", "latitude", 
                        "longitude", 'type', "quality", "error_semi_major_axis", 
                        "error_semi_minor_axis", "error_ellipse_orientation"))
      
      # # Next, join the data; first we pare down the data frames a bit
      tbl_lsub <- dplyr::select(locs, c("deployid", "datetime", "latitude",
                                        "longitude", "type", "quality", "error_semi_major_axis",
                                        "error_semi_minor_axis", "error_ellipse_orientation"))
      
      daflocs <- dplyr::bind_rows(deploy_dat, tbl_lsub) %>% 
        arrange(deployid, datetime)
      daflocs$deployid <- stringr::str_remove(daflocs$deployid, "_DUML")
      
      
      # Goniometer Data
      # Add a comment that gonio_data is a target
      # is there a style guide for this
      gonio <- gonio_data %>%
        dplyr::mutate(quality = 'User',
               type = 'User') %>% 
        dplyr::distinct() %>% 
        janitor::clean_names() %>%
        dplyr::rename(datetime = 'date', 
                      latitude = 'lat', longitude = 'lon') %>% 
        dplyr::arrange(deployid) 
      
      gonio <- gonio %>%
        dplyr::filter(!gps_gonio_time_difference_s > 60, strength_db <= 0)
      gonio$gap <- c(0, diff(gonio$datetime) > 120) 
      gonio$group <- cumsum(gonio$gap) + 1
      
      gonio_out <- gonio %>% 
        dplyr::group_by(deployid, group) %>% 
        dplyr::filter(strength_db == max(strength_db))
      
      gonio_out <- gonio_out %>% 
        dplyr::group_by(deployid, group) %>% 
        dplyr::slice_head() %>%
        ungroup()
      
      # Tncorporating the curve from Anna Borroni
      # > als <- coefficients(lm(log10(Gonio_dB_distance_curve$distance) ~Gonio_dB_distance_curve$dB ))
      # > als
      # (Intercept) Gonio_dB_distance_curve$dB 
      # 1.06103672                -0.02690245 
      # 10^(als[1] + (als[2] * -120))
      anna_params <- c(1.06103672, -0.02690245)
      gonio_out <- gonio_out %>% 
        dplyr::mutate(distance = 10^(anna_params[1] + (anna_params[2] * strength_db))) %>%
        dplyr::select(deployid, datetime, latitude, longitude, type, quality, platform, distance) %>%
        dplyr::mutate(quality = rep("User", nrow(gonio_out)),
               source = rep("gonio", nrow(gonio_out))) 
      
      gonio_out <- gonio_out %>% 
        dplyr::mutate(error_semi_major_axis = distance,
               error_semi_minor_axis = distance,
               error_ellipse_orientation = 0) %>% 
        dplyr::select(deployid, datetime, latitude, longitude, type, quality, 
               error_semi_major_axis, error_semi_minor_axis, error_ellipse_orientation)
      
      daflocs2 <- bind_rows(daflocs, gonio_out) %>% 
        dplyr::group_by(deployid) %>% 
        dplyr::arrange(datetime)
      
      # Remove the times before one day ahead of the deployment
      daflocs2 <- dplyr::inner_join(daflocs2, deploy_dat, 
                                    by = "deployid", 
                                    suffix = c("_location", "_deploy"))
      
      daflocs2 <- daflocs2 %>%
        dplyr::filter(datetime_location >= datetime_deploy - lubridate::days(1))
      
      daflocs2 <- daflocs2 %>%
        dplyr::rename_with(~ str_remove(., "_location"), 
                    ends_with("_location")) %>%
        dplyr::select(-ends_with("_deploy"))
      
      # Remove duplicated times
      locs <- daflocs2 %>%
        dplyr::filter(deployid %in% exposed_ans) %>% 
        dplyr::group_by(deployid) %>%
        dplyr::arrange(datetime, error_semi_major_axis) %>%
        tidyr::drop_na(error_semi_major_axis) %>% 
        tidyr::drop_na(longitude) %>% 
        dplyr::mutate(
          rank = 1L,
          rank = case_when(duplicated(datetime, fromLast = FALSE) ~
                             lag(rank) + 1L, TRUE ~ rank)) %>%
        dplyr::filter(rank == 1) %>%
        dplyr::arrange(deployid, datetime) %>%
        readr::write_csv(file = here::here('data/proc_data',
                                           'ss_cleaned_locs.csv'))
      
      return(locs)
      
      
    }
  )
)