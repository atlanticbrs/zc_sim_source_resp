register_target(
  
  # Here we do the sampling of the RLs
  # I need to read in the imputed points: x,y,z,t: pts_mult_imp_rejsamp
  # the RL manifest
  
  tar_target(
    
    name = calc_rl_rj,
    command = {
      
      # Container to hold output
      rls_rj <- vector(mode = 'list', length = length(pts_mult_imp_rejsamp))
      # loop over one animal/cee combination at a time
      for(i in 1:length(pts_mult_imp_rejsamp)){
        
        my_cee_id <- pts_mult_imp_rejsamp[[i]][["cee_id"]]
        animal <- pts_mult_imp_rejsamp[[i]][["imputed_points"]] %>% 
          sf::st_as_sf() %>% 
          mutate(cee_id = my_cee_id,
                 d2ship_km = NA,
                 bearing = NA,
                 rl_lci = NA,
                 rl_med = NA,
                 rl_uci = NA,
                 all_rl_vals = NA,
                 tobl = NA)
        
        sub_manifest <- data.frame(rl_dat[rl_dat$cee_id == my_cee_id, ])
        
        # Within the animal/cee combo, loop over the 
        # cee times to extract RL values
        for(j in 1:nrow(sub_manifest)){
          
          # Path to the RL netCDF File
          zc_foot <- nc_open(sub_manifest[j, 'paths'])
          
          # Vector of times over which we'll extract RLs
          tvec <- sub_manifest$rl_time[j]
          
          # lat/lon of the boat
          boat_loc_ll <- sub_manifest[j, c('lon', 'lat')]
          boat_loc_ea <- boat_loc_ll %>% 
            sf::st_as_sf(coords = c("lon", "lat")) %>%
            sf::st_set_crs(4326) %>% 
            sf::st_transform(brs_crs)
          
          # Calculate Distance to the Ship for each Imputed Point at the right time
          idx <- which(animal$datetime == tvec)
          d2ship_m <- sf::st_distance(animal[idx, ], boat_loc_ea)
          animal[idx, 'd2ship_km'] <- d2ship_m / 1000
          
          # Calculate Bearing to Ship for each Imputed Point
          my_df_ll <- animal %>% 
            sf::st_transform(crs = 4326) %>% 
            sf::st_coordinates()
          
          # Use the actual bearing to get a high index and a low index so I can do 
          # interpolation between radials
          # Have to error trap for range between 0 and 1, because 
          # ceiling yields 0, which yields an error
          my_bearing <- geosphere::bearingRhumb(boat_loc_ll, my_df_ll[idx, ])
          animal[idx, 'bearing_low'] <- floor(my_bearing) 
          animal[idx, 'bearing_high'] <- ceiling(my_bearing)
          
          # get the RL Info
          zc_lat <- ncvar_get(zc_foot, 'latitude') # dimensions of range
          zc_lon <- ncvar_get(zc_foot, 'longitude') # dimensions as above for lat
          zc_rl <- ncvar_get(zc_foot, 'RL') # dimensions are[range, radial, depth]
          zc_bearing <- ncvar_get(zc_foot, 'bearing') # id of each of the radials
          zc_range <- ncvar_get(zc_foot, 'range') # units are kilometers
          zc_depth <- ncvar_get(zc_foot, 'depth') # Just the bins, i.e. 0-10, 11-20, etc.
          
          ### Extract the RLs
          # idx is the vector of points for this ith 5-minute period
          # in the cee
          for(k in idx){
            
            # The RL data in zc_rl is an array with dimensions:
            # 1. depth
            # 2. radial
            # 3. range from ship
            # We need to find these indices to extract the RL data
            # range_bin gets distance from ship
            # blowidx and bhighidx gets used to locate the right Radial slices
            # if blowidx == 0, i.e. True North, then set it to 360
            # dminidx & dmaxidx get the depth bins in the water columm
            range_bin <- findInterval(animal$d2ship_km[k], zc_range)
            blowidx <- animal$bearing_low[k]
            if(blowidx == 0) blowidx <- 360
            bhighidx <- animal$bearing_high[k]
            dminidx <- findInterval(animal$depth_min[k], zc_depth)
            dmaxidx <- findInterval(animal$depth_max[k], zc_depth)
            
            # Here we get the range of index values for this animal's depth bin range
            # if NAs, take the whole water column
            if(any(is.na(c(dminidx, dmaxidx))) | (dminidx == 0 & dmaxidx == 0)){
              watercolidx <- 1:length(zc_depth)
            } else {
              watercolidx <- dminidx:dmaxidx  
            }
            
            # Extract & summarize
            # sound is a 3 column data frame, whose # of rows corresponds to the 
            # depth min to depth max range for a single point
            sound <- data.frame(rl_low = zc_rl[range_bin, blowidx, watercolidx],
                                rl_high = zc_rl[range_bin, bhighidx, watercolidx]) 
            rl_mean <- rowMeans(sound, na.rm = TRUE)
            
            animal$rl_lci[k] <- quantile(rl_mean, probs = c(0.025, .5, 0.975), na.rm = TRUE)[1]
            animal$rl_med[k] <- quantile(rl_mean, probs = c(0.025, .5, 0.975), na.rm = TRUE)[2]
            animal$rl_uci[k] <- quantile(rl_mean, probs = c(0.025, .5, 0.975), na.rm = TRUE)[3]
            animal$all_rl_vals[k] = list(rl_mean)
            
          } # end k loop over k - the index of imputed points 
          
        } # end j loop over the cee times in sub_manifest
        
        # Make An exploratory plot
        # p <- ggplot2::ggplot(animal, aes(x = rl_med, y = factor(datetime), fill = after_stat(x)))+
        #   ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
        #   theme_bw()+
        #   labs(title = "Estimated RLs at 5' Intervals within the CEE",
        #        subtitle = paste0(unique(animal$deployid), ", CEE ID: ", my_cee_id),
        #        y = '', x = 'dB RMS',
        #        fill = 'dB RMS')
        #   scale_fill_viridis_c()
        #   
        #   ggsave(p, filename = here::here("results/plots", 
        #                                paste0(unique(animal$deployid), "_", my_cee_id, "_rl-plot.png")),
        #          device = 'png',
        #          dpi = 'retina',
        #          width = 10, height = 6.1, units = 'in')  
        
        # Write Animals to Disk - animal is both an sf object and data frame
        saveRDS(animal, file = here::here('data/proc_data', 
                                          paste0(unique(animal$deployid), "_", my_cee_id, '_rl_rj.rds')))
        
        # Store items in List 
        rls_rj[[i]] <- list(deployid = unique(animal$deployid),
                            cee_id = my_cee_id,
                            pts_rl = animal)
        
        
      } # end i loop over the list of imputed points for each animal/cee combination
      
      rls_rj
      
      })
             

)