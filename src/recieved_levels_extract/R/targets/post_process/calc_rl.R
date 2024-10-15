register_target(
  
  # Here we do the sampling of the RLs
  # I need to read in the imputed points: x,y,z,t: pts_mult_imp
  # the RL manifest
  
  tar_target(
    
    name = calc_rl,
    command = {
      
      # Container to hold output
      rls <- vector(mode = 'list', length = length(pts_mult_imp))
      # loop over one animal/cee combination at a time
      for(i in 1:length(pts_mult_imp)){
        
        my_cee_id <- pts_mult_imp[[i]][["cee_id"]]
        outname <- pts_mult_imp[[i]][["deployid"]]
        my_cee_times <- data.frame(pts_mult_imp[[i]]['pred_times'])
        series_data <- pts_mult_imp[[i]][["series_stats"]]
        series_cee_gaps <- pts_mult_imp[[i]][["series_gaps"]]
        animal <- pts_mult_imp[[i]][["predicted_pts"]] %>% 
          sf::st_as_sf() %>% 
          mutate(cee_id = my_cee_id,
                 d2ship_km = NA,
                 bearing = NA,
                 rl_lci = NA,
                 rl_med = NA,
                 rl_uci = NA,
                 all_rl_vals = NA,
                 tobl = NA)
        
        # read in the data from water_col_idx
        sub_wcol <- an_series_df %>% 
          filter(deployid == outname,
                 cee_id == my_cee_id)
        
        print(paste(outname, my_cee_id, sep = " - "))
        print(paste("Series data?", series_data, sep = " - "))
        print(paste("Gaps in Series Data During CEE?", series_cee_gaps, sep = " - "))
        
        # Test if animal was exposed to sound, i.e. not a control
        exposed <- my_cee_id %in% rl_dat$cee_id
        
        if(exposed){
          
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
              range_bin <- findInterval(animal$d2ship_km[k], zc_range)
              blowidx <- animal$bearing_low[k]
              if(blowidx == 0) blowidx <- 360
              bhighidx <- animal$bearing_high[k]
              
              # Here we take extract the indices for the water column
              if(series_data & outname != "ZcTag073_DUML" & is.na(series_cee_gaps)){
                
                # assemble indices
                dminidx <- findInterval(sub_wcol$depth_min[j], zc_depth)
                dmaxidx <- findInterval(sub_wcol$depth_max[j], zc_depth)
                
                # Here we get the range of index values for this animal's depth bin range
                # if NAs, take the whole water column
                if(any(is.na(c(dminidx, dmaxidx))) | (dminidx == 0 & dmaxidx == 0) ){
                  watercolidx <- 1:length(zc_depth)
                } else{
                  watercolidx <- dminidx:dmaxidx  
                }
                
              } else { # if exposed, but no series data or gaps
                watercolidx <- 1:length(zc_depth)
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
          
        } else { # if not exposed, just calculate distance to ship
          
          for(j in 1:nrow(my_cee_times)){
            
            # Vector of times over which we'll extract RLs
            tvec <- my_cee_times[j, ]
            
            # lat/lon of the boat
            boat_loc_ll <- cee_dat %>% 
              filter(cee_id == my_cee_id) %>% 
              slice_head() %>% 
              dplyr::select("lon_start", "lat_start")
              
            boat_loc_ea <- boat_loc_ll %>% 
              sf::st_as_sf(coords = c("lon_start", "lat_start")) %>%
              sf::st_set_crs(4326) %>% 
              sf::st_transform(brs_crs)
            
            # Calculate Distance to the Ship for each Imputed Point at the right time
            idx <- which(animal$datetime == tvec)
            d2ship_m <- sf::st_distance(animal[idx, ], boat_loc_ea)
            animal[idx, 'd2ship_km'] <- d2ship_m / 1000
            
            
          } # end j loop over the cee times 
          
        } # end distance calculations
        
        # Write Animals to Disk - animal is both an sf object and data frame
        saveRDS(animal, file = here::here('data/proc_data', 
                                          paste0(unique(outname), "_", my_cee_id, '_rl.rds')))
        
        # Store items in List 
        rls[[i]] <- list(deployid = unique(outname),
                            cee_id = my_cee_id,
                            pts_rl = animal)
        
        
      } # end i loop over the list of imputed points for each animal/cee combination
      
      rls
      
      })
             

)