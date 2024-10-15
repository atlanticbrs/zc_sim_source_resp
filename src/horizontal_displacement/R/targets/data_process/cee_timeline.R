register_target(
  # Make a line plot of the deployment
  tar_target(
    name = cee_timeline,
    command = {
      
      deploy_timeline <- cee_dat %>% 
        sf::st_drop_geometry() %>%
        dplyr::select("deployid", "cee_id", "loc_st", "loc_en") %>% 
        mutate(deployid = gsub('Zc', 'S', deployid),
               deploy_start_time = as.POSIXct(loc_st, origin = "1970-01-01", tz = 'UTC'),
               deploy_end_time = as.POSIXct(loc_en, origin = "1970-01-01", tz = 'UTC'),
               year = lubridate::year(deploy_start_time),
               winyear_start = format(deploy_start_time, "%m-%d %H:%M:%S"),
               winyear_start = as.POSIXct(winyear_start, format = "%m-%d %H:%M:%S", tz = "UTC"),
               winyear_end = format(deploy_end_time, "%m-%d %H:%M:%S"),
               winyear_end = as.POSIXct(winyear_end, format = "%m-%d %H:%M:%S", tz = "UTC")) %>% 
        dplyr::select(-c(deploy_start_time, deploy_end_time, "loc_st", "loc_en")) %>% 
      dplyr::group_by(deployid) %>% 
        slice_head() 
        
      
      cee_times <- cee_dat %>% 
        sf::st_drop_geometry() %>%
        dplyr::select("deployid", "cee_id", "cee_type", "cee_st") %>% 
        dplyr::mutate(deployid = gsub('Zc', 'S', deployid),
               cee_start_time = as.POSIXct(cee_st, origin = "1970-01-01", tz = 'UTC'),
               winyear_cee_start = format(cee_start_time, "%m-%d %H:%M:%S"),
               winyear_cee_start = as.POSIXct(winyear_cee_start, format = "%m-%d %H:%M:%S", tz = "UTC"),
               year = lubridate::year(cee_start_time)) %>% 
        dplyr::select(-cee_st)
      
      pbase <- ggplot(deploy_timeline)+
        geom_segment(aes(x = winyear_start, xend = winyear_end, y = deployid), colour = '#525252')+
        theme_bw()+
        scale_x_datetime(date_labels = "%b", date_breaks = "1 month") +
        # facet_grid(vars(year), scales = 'free_y', space = 'free_y')
        facet_wrap(vars(year),ncol = 3, scales = 'free_y')+
        labs(x = 'Date', y = 'Deploy ID')
      
      p_points <- pbase +
        geom_point(data = cee_times, aes(x = winyear_cee_start, y = deployid, color = cee_type))+
        facet_wrap(vars(year),ncol = 3, scales = 'free')+
        scale_color_brewer(palette = 'Dark2', type = 'qual')+
        labs(color = 'CEE Type')
      
      # p_shift <- shift_legend2(p_points)
      
      f = file.path(here::here('results'))
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      ggsave(p_points, filename = file.path(f, paste(Sys.Date(), '_', 'cee_extents.png', sep ='')),
             device = 'png',
             dpi = 'retina',
             width = 10, height = 6.1, units = 'in') 

    }
  )
)