register_target(
  # Make the CEE data spatial
  tar_target(
    name = cee_sf,
    command = {
      
      cee_sf <- cee_dat %>%
        dplyr::group_by(cee_id) %>% 
        dplyr::slice_head() %>% 
        sf::st_as_sf(coords = c("lon_start","lat_start")) %>%
        sf::st_set_crs("EPSG:4326")
      
      cee_sf <- cee_sf %>% 
        sf::st_transform(brs_crs)
      
      return(cee_sf)
    }
  )
)