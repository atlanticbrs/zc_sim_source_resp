register_target(
  
  # take the locations (that are sf objects) and convert
  # them to a trip object in order to run the sda filter
  # from the trip library
  
  tar_target(
    name = locs_filter,
    command = {
      
      crs <- sf::st_crs(locs_sf)
      dat <- locs_sf %>%
        sf::st_transform(4326) %>%
        ungroup() %>%
        arrange(deployid, datetime)
      dat_tr <- trip::trip(dat, c("datetime", "deployid"), correct_all = FALSE)
      
      # Employ Speed Filter - using trip library
      # 2.7778 m/s for Ziphius
      # 10 km/hour for Ziphius and 15 km/hr for pilots 
      suppressWarnings(
        keep <- sda(
          dat_tr,
          smax = 10 #km/hour or 8.33m/sec * 3.6
        )
      )
      
      tracks_filt <- dat %>%
        mutate(sda_keep = keep) %>%
        filter(sda_keep) %>%
        dplyr::select(-c(sda_keep, rank)) %>%
        st_transform(crs)
      
      return(tracks_filt)
      
    }
  )
  
  
)