register_target(
  # 
  tar_target(
    name = prep_rasters,
    command = {
      
      # read in the Bathy raster
      bathy_file <- here::here("data/Extract_topo30Rectangle_Proj.tif")
      bathy <- raster::raster(bathy_file)
      bathy <- raster::disaggregate(bathy, fact = 2, method = "bilinear") # results in 500m grid
      ch_stack <- raster::stack(bathy)
      names(ch_stack) <- 'bathy'
      
      # Make a shelf-break/off-shelf raster for location covariate
      m1 <- c(-6000, -200, 1,  -200, 2000, 0)
      rclmat <- matrix(m1, ncol=3, byrow=TRUE)
      shelf_brk <- raster::reclassify(bathy, rclmat)
      names(shelf_brk) <- 'shelf_brk'
      ch_stack <- raster::addLayer(ch_stack, shelf_brk)
      
      # Create cee-specific distance to ship rasters
      for(i in 1:nrow(cee_sf)){
        
        cee_sub <- cee_sf[i, ]
        
        d2ship <- bathy
        d2ship <- raster::distanceFromPoints(d2ship, st_coordinates(cee_sub))
        names(d2ship) <- paste0("d2ship_cee_", cee_sub$cee_id)
        
        ch_stack <- raster::addLayer(ch_stack, d2ship)
        
      }
      
      ch_stack
    }
  )
)