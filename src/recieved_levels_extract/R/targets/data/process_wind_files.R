register_target(
  
  # Here we:
  # 1. read in the wind files one at a time, 
  # 2. intersect them with the cee information to get the time slices
  # 3. pull out the right wind raster
  # 4. save it to disk
  
  tar_target(
    
    name = process_wind_files,
    command = {
      
      start_wind_times <- c(strptime("2017-05-01 00:00:00", 
                                     format = "%Y-%m-%d %H:%M:%S", tz = 'UTC'),
                            strptime("2018-05-01 00:00:00", 
                                     format = "%Y-%m-%d %H:%M:%S", tz = 'UTC'),
                            strptime("2019-05-01 00:00:00", 
                                     format = "%Y-%m-%d %H:%M:%S", tz = 'UTC'),
                            strptime("2020-05-01 00:00:00", 
                                     format = "%Y-%m-%d %H:%M:%S", tz = 'UTC'))
      
      # Container to hold the information about wind data for each cee
      # Note that the last two in 2022 are controls - hence
      # subtracting off 2
      wind_list <- vector(mode = 'list', length = length(cee_times) - 2)
      
      for(i in 1:length(wind_list)){
        
        cee_time <- min(cee_times[[i]]$times)
        outname <- cee_times[[i]][["cee_id"]]
        yidx <- lubridate::year(cee_time)
        
        fidx <- case_when(
          yidx == 2017 ~ 1,
          yidx == 2018 ~ 2,
          yidx == 2019 ~ 3,
          yidx == 2020 ~ 4
        )
        
        # generic code
        wind_file = wind_files_directory[fidx]
        start_wind_time <- start_wind_times[fidx]
        
        w_foot <- nc_open(wind_file)
        w_date_idx <- ncvar_get(w_foot, 'Date') 
        w_lon <- ncvar_get(w_foot, 'Longitude')[, 1]
        w_lat <- ncvar_get(w_foot, 'Latitude')[1, ]
        
        w_date <- seq.POSIXt(from = start_wind_time, 
                             by = '3 hours',
                             length.out = length(w_date_idx))
        
        cee_idx <- min(which(w_date > cee_time))
        
        # pull out the nearest good day with TOBL values
        # for three different cees in 2018
        if(outname == "18_02") cee_idx <- 249
        if(outname == "18_03") cee_idx <- 249
        if(outname == "18_08") cee_idx <- 985
        
        # Raster Processing: spatial extent, temporal extent,
        # projection, and removing NAs
        wind_nc <- raster::brick(wind_file)
        raster::extent(wind_nc) <- raster::extent(min(w_lon), max(w_lon), 
                                  min(w_lat), max(w_lat))
        wind_nc_cee <- wind_nc[[cee_idx]]
        wind_nc_cee <- raster::clamp(wind_nc_cee, lower = 0)
        wind_nc_cee_prj <- raster::projectRaster(wind_nc_cee, crs = brs_crs)
        
        # Plot for examining
        # png(here::here("results/plots/wind-files", 
        #                paste0(cee_time, "_wind-plot.png")))
        # raster::plot(wind_nc_cee, main = cee_time)
        # dev.off()
        
        # save to disk
        outfilename <- here::here('data/proc_data', 
                                  paste0('wind_raster_ceeid_', outname, ".grd"))
        raster::writeRaster(wind_nc_cee_prj, 
                filename = outfilename, overwrite=TRUE)
        
        wind_list[[i]] <- data.frame(cee_id = outname,
                               cee_time = cee_time,
                               wrast_path = outfilename)
        
        
      }
      
      wind_df <- do.call('rbind', wind_list)
      # data.frame(wind_df)
      
      })
             

)