library(tidyverse)
library(here)
library(readr)

# All Tags
tags_all <- read_csv("data/tagstreams_by_cee_details_filt.csv",
                     show_col_types = FALSE) %>% 
  filter(cee_type %in% c("simulated mfas", "control"),
         is.na(ser_gaps)) %>% 
  filter(ser == TRUE) %>% 
  filter(deployid != "ZcTag079_DUML") %>% 
  mutate(baseline_end = as.POSIXct.numeric(cee_st, origin = "1970-01-01", tz = 'UTC'),
         cee_start = as.POSIXct.numeric(cee_st, origin = "1970-01-01", tz = 'UTC'),
         series_start = as.POSIXct.numeric(ser_st, origin = "1970-01-01", tz = 'UTC'),
         series_end = as.POSIXct.numeric(ser_en, origin = "1970-01-01", tz = 'UTC'),
         sex = "U",
         deployid = gsub(pattern = '_DUML', replacement = '', x = deployid))  %>% 
  select(deployid, sex, baseline_end, cee_start, cee_id, series_start, series_end) 

# this is a good time before a gap for ZcTag073
# this is a tag that reset, so has a long series record, but with big gaps
# that mess up the pre-post indexing
zctag073_end <- as.POSIXct("18-Aug-2018 23:55:00", format = "%d-%b-%Y %H:%M:%S", tz = 'UTC')
tags_all[tags_all$deployid == 'ZcTag073', 'series_end'] <- zctag073_end

# set up directory for dive data
dive_dir <- "../../01_shared_data_products/filter_sattag/series/"
my_files <- data.frame(my_full_files = dir(dive_dir,
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE,
                                           full.names = TRUE),
                       my_file_names = dir(dive_dir,
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE))
my_files$deployid <- stringr::str_extract(my_files$my_file_names, 'ZcTag[0-9]+')

tags_all$ser_length <- NA
tags_all$cee_start_idx<- NA
for(i in 1:nrow(tags_all)){
  
  fidx <- which(tags_all$deployid[i] == my_files$deployid)
  animal <- read_csv(my_files$my_full_files[fidx], show_col_types = FALSE)
  an_date_time <- animal %>% 
    pluck('Date')
  an_date_time <- as.POSIXct.numeric(an_date_time, origin = "1970-01-01", tz = 'UTC')
  tags_all$ser_length[i] <- nrow(animal)
  
  idx <- min(which(tags_all$baseline_end[i] <= an_date_time))
  tags_all$cee_start_idx[i] <- idx
  
  if(tags_all$deployid[i] == 'ZcTag073'){
    tseq <- seq(tags_all$series_start[i], tags_all$series_end[i], by = '5 mins')
    tags_all$ser_length[i] <- length(tseq)
  }
  
}

write_csv(tags_all, file = here::here('data/tag_info_all-tags_manifest.csv'))

# Assemble the baseline segments
tags_mult <- tags_all %>% 
  mutate(baseline_start_idx = 1)

# Write out CEE Start indices. Last CEE in the record
# must be at least 6 hours before the end of the data
# else there's not enough data to run the pre-post
tags_mult_cee_out <- tags_mult %>% 
  filter(ser_length - cee_start_idx > 6 * 12,
         cee_start_idx > 6 * 12) %>% 
  select(deployid, cee_start_idx)#, cee_id
write_csv(tags_mult_cee_out, file = here::here('data/sattag/tag_info_mult-segments_cee.csv'))

# Calculate and Write-out the Baseline Indices
uniq_ans <- unique(tags_mult_cee_out$deployid)
tags_mult_bline_out <- numeric()
for(my_an in uniq_ans){
  df <- tags_mult %>% 
    filter(deployid == my_an)%>% 
    mutate(baseline_end_idx = cee_start_idx)
  
  if(nrow(df) == 1){
    
    last_row <- df %>% 
      mutate(baseline_start_idx = cee_start_idx + 24*12,
             baseline_end_idx = ser_length)
    df <- bind_rows(df, last_row)
    
  } else {
    
    for(i in 2:nrow(df)){
      
      # Add a 24-hour period after a CEE
      df$baseline_start_idx[i] <- df$cee_start_idx[i-1] + 24*12
      
    } # end inner loop over animal
    
    last_row <- df %>% 
      slice_tail %>% 
      mutate(baseline_start_idx = baseline_end_idx + 24*12,
             baseline_end_idx = ser_length)
    df <- bind_rows(df, last_row)
    
  }
  
  tags_mult_bline_out <- rbind(tags_mult_bline_out, df)
  
}

tags_mult_bline_out <- tags_mult_bline_out %>% 
  filter(ser_length - baseline_start_idx > 0,
         ser_length - (cee_start_idx + (6 * 12)) > 0,
         cee_start_idx > 6 * 12)

# Assemble a full record to have on hand for reference
tags_mult_bline_metadata <- tags_mult_bline_out 

tags_mult_bline_out <- tags_mult_bline_out %>% 
  select(deployid, baseline_start_idx, baseline_end_idx)

write_csv(tags_mult_bline_out, file = here::here('data/tag_info_mult-segments_baseline.csv'))
write_csv(tags_mult_bline_metadata, file = here::here('data/tag_info_mult-segments_baseline_metadata.csv'))


