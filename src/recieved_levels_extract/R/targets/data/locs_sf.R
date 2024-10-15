register_target(
  # Convert the tracks to an sf object
  
  tar_target(
    name = locs_sf,
    command = {
      
      tracks_sf <- locs %>%
        sf::st_as_sf(coords = c("longitude","latitude")) %>%
        sf::st_set_crs("EPSG:4326")
      
      tracks_sf <- tracks_sf %>% 
        sf::st_transform(brs_crs)
      
      return(tracks_sf)
    }
  )
)  