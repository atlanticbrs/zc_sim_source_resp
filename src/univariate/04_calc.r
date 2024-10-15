###
# calc.r
# calculate metrics for all applicable tags 



### constants
EXPOSURE_WINDOW <- 24 * 60 * 60 # 1 day (in seconds)
DEEP_DIVE_THRESHOLD <- 800 # meters



### load data
load("02_intermediate_outputs/setup_define_out.rdata")



### calculate metrics
# (S) = series only

# 1. duration of exposure dive 
# 2. (S) ascent duration of exposure dive
# 3. depth of exposure dive
# 4. inter deep dive interval (IDDI) following exposure dive
# 5. inter deep dive interval (IDDI) prior to the exposure dive
# 6. (S) n dives in the IDDI
# 7. (S) mean dive duration during the IDDI
# 8. (S) mean max depth of dives during the IDDI
# 9. (S) max depth during the IDDI

## other details:
# cee commense + 24 hours is the non baseline section
# everything else is baseline

# list to hold the data
plot_inputs <- list()
count <- 1 # counter to keep track of place in list

# megaloop
for(i in 1:nrow(scaled)) {
  ceest <- scaled$st_datenum[i]
  ceeen <- scaled$en_datenum[i]
  
  # behavior
  for(b in 1:length(behavior_bout_summary)) {
    bb_cur <- behavior_bout_summary[[b]]
    
    # get index of exposure dive
    exp_ind <- which(ceest >= bb_cur$dive_st & ceest <= bb_cur$dive_en)
    
    # keep track if exposure happened during a dive or a surf
    surf_ind <- FALSE
    
    # try exposure IDDI
    if(length(exp_ind) == 0) {
      exp_ind <- which(ceest >= bb_cur$surf_st & ceest <= bb_cur$surf_en)
      surf_ind <- TRUE
    }

    # dives are exclusive to time so there should never be more than one match
    if(length(exp_ind) > 0) {
      if(length(exp_ind) > 1) {
        error(paste(i, b, "beahvior. more than one matching exposure dive"))
      }
      
      # generate baseline
      after_exposure_baseline_st <- ceest + EXPOSURE_WINDOW # seconds
      baseline <- bb_cur$dive_en < ceest | bb_cur$dive_st >= after_exposure_baseline_st
      
      # select deep or shallow dives
      selected_baseline <- baseline & bb_cur$depth > DEEP_DIVE_THRESHOLD
      label_selected <- "exp. dive > 800 meters"
      
      # if the exposure dive is shallow only compare to shallow
      if(bb_cur$depth[exp_ind] <= DEEP_DIVE_THRESHOLD | surf_ind) {
        selected_baseline <- baseline & bb_cur$depth <= DEEP_DIVE_THRESHOLD
        label_selected <- "exp. dive \u2264 800 meters"
      }
      
      # specify baseline if exp is in the surface
      if(surf_ind) {
        selected_baseline <- baseline
      }
      
      # gather info
      info <- list(
        ceeid = scaled$cee_id[i],
        ceetype = scaled$cee_type[i],
        ceest = ceest,
        ceeen = ceeen,
        dat_raw = beh[which(sapply(behavior_bout_summary, function(x) x$DeployID[1]) == bb_cur$DeployID[1])],
        dat_pro = bb_cur, 
        deployid = bb_cur$DeployID[1],
        label_selected = label_selected,
        settings = "behavior"
      )
      
      if(!surf_ind) {
      # output dive duration
      plot_inputs[[count]] <- list(
        obs = bb_cur$dive_duration[exp_ind] / 60, # in min observed exposure dive dur
        bas = bb_cur$dive_duration[selected_baseline] / 60, #in min baseline dive dur
        info = info,
        what = "dive_duration"
      )
      count <- count + 1
      
      
      # output dive max depth
      plot_inputs[[count]] <- list(
        obs = bb_cur$depth[exp_ind], # meters
        bas = bb_cur$depth[selected_baseline], # meters
        info = info,
        what = 'dive_maxdepth'
      )
      count <- count + 1
      }
      
      # output IDDI if there is one
      if(!is.na(bb_cur$surf_duration[exp_ind])) {
        
        mod_baseline <- selected_baseline
        
        # modify baseline so it doesn't include the previous IDDI either
        # since this IDDI can be interrupted by the exposure
        # unless the exposure happens in an IDDI
        if(!surf_ind) {
        if(!is.na(bb_cur$surf_duration[exp_ind-1])) {
          mod_baseline <- (1:length(selected_baseline) != exp_ind - 1) & selected_baseline
        }}
        
        plot_inputs[[count]] <- list(
          obs = bb_cur$surf_duration[exp_ind] / 60, # in minutes observed IDDI
          bas = bb_cur$surf_duration[mod_baseline] / 60, # in minutes baseline IDDI
          info = info,
          what = 'IDDI'
        )
        count <- count + 1
      }
      
      # output the prior IDDI if there is one
      if(!is.na(bb_cur$surf_duration[exp_ind-1])) {
        plot_inputs[[count]] <- list(
          obs = bb_cur$surf_duration[exp_ind - 1] / 60, # in minutes observed prior IDDI
          bas = bb_cur$surf_duration[selected_baseline] / 60, # in minutes baseline IDDI
          info = info,
          what = 'prior_IDDI'
        )
        
        count <- count + 1
      }
    }
  }
  
  # series dive summary
  # exp dive duration
  # mas depth
  # ascent duration
  for(s in 1:length(series_dive_summary)) {
    sd_cur <- series_dive_summary[[s]]
    
    # get index of exposure dive (if there is one)
    exp_ind <- which(ceest >= sd_cur$Start & ceest <= sd_cur$End)
    
    # dives are exclusive to time so there should never be more than one match
    if(length(exp_ind > 0)) {
      if(length(exp_ind) > 1) {
        error(paste(i, s, "series. more than one matching exposure dive"))
      }
      
      # generate baseline
      after_exposure_baseline_st <- ceest + EXPOSURE_WINDOW # seconds
      baseline <- sd_cur$End < ceest | sd_cur$Start >= after_exposure_baseline_st
      
      # select deep or shallow dives
      selected_baseline <- baseline & sd_cur$maxdepth > DEEP_DIVE_THRESHOLD
      label_selected <- "exp. dive > 800 meters"
      
      # if the exposure dive is shallow only compare to shallow
      if(sd_cur$maxdepth[exp_ind] <= DEEP_DIVE_THRESHOLD) {
        selected_baseline <- baseline & sd_cur$maxdepth <= DEEP_DIVE_THRESHOLD
        label_selected <- "exp. dive \u2264 800 meters"
      }
      
      info <- list(
        ceeid = scaled$cee_id[i],
        ceetype = scaled$cee_type[i],
        ceest = ceest,
        ceeen = ceeen,
        dat_raw = ser[which(sapply(series_bout_summary, function(x) x$DeployID[1]) == sd_cur$DeployID[1])],
        dat_pro = sd_cur, 
        deployid = sd_cur$DeployID[1],
        label_selected = label_selected,
        settings = "series"
      )
      
      # output dive_duration
      plot_inputs[[count]] <- list(
        obs = sd_cur$nsamp[exp_ind] * 5, # in min observed exposure dive dur
        bas = sd_cur$nsamp[selected_baseline] * 5, # in min baseline dive dur
        info = info,
        what = "dive_duration"
      )
      count <- count + 1
      
      # output dive_maxdepth
      plot_inputs[[count]] <- list(
        obs = sd_cur$maxdepth[exp_ind], # meters
        bas = sd_cur$maxdepth[selected_baseline], # meters
        info = info,
        what = "dive_maxdepth"
      )
      count <- count + 1
      
      # output ascent duration
      plot_inputs[[count]] <- list(
        obs = sd_cur$nascent[exp_ind] * 5, # in minutes exposure
        bas = sd_cur$nascent[selected_baseline] * 5, # in minutes baseline
        info = info,
        what = "dive_ascent_duration"
      )
      count <- count + 1
    }
  }
  
  # series bout summary
  #  n dives in the IDDI
  #  mean dive duration during the IDDI
  #  mean max depth of dives during the IDDI
  #  max depth during the IDDI
  for(sb in 1:length(series_bout_summary)) {
    sb_cur <- series_bout_summary[[sb]]
    
    # get index of exposure dive (if there is one)
    exp_ind <- which(ceest >= sb_cur$Start & ceest <= sb_cur$End)
    
    # dives are exclusive to time so there should never be more than one match
    if(length(exp_ind > 0)) {
      if(length(exp_ind) > 1) {
        error(paste(i, sb, "series bout. more than one matching exposure dive"))
      }
      
      # generate baseline
      after_exposure_baseline_st <- ceest + EXPOSURE_WINDOW # seconds
      baseline <- sb_cur$End < ceest | sb_cur$Start >= after_exposure_baseline_st
      
      # select deep or shallow dives
      label_selected <- "exp. dive > 800 meters"
      
      # if the exposure dive is shallow indicate it
      sd_cur <- series_dive_summary[[sb]]
      sd_exp_ind <- which(ceest >= sd_cur$Start & ceest <= sd_cur$End)
      if(sd_cur$maxdepth[sd_exp_ind] <= DEEP_DIVE_THRESHOLD) {
        label_selected <- "exp. dive \u2264 800 meters"
      }
      
      info <- list(
        ceeid = scaled$cee_id[i],
        ceetype = scaled$cee_type[i],
        ceest = ceest,
        ceeen = ceeen,
        dat_raw = ser[which(sapply(series_bout_summary, function(x) x$DeployID[1]) == sd_cur$DeployID[1])],
        dat_pro = sd_cur, 
        deployid = sb_cur$DeployID[1],
        label_selected = label_selected,
        settings = "series"
      )
      
      # modified baseline for IDDI just in case prior IDDI is interrupted by exp
      modified_baseline <- baseline
      if(exp_ind - 1 > 0) {
        modified_baseline <- (1:length(baseline) != exp_ind - 1) & baseline
      }
      
      
      # output IDDI
      iddi_dur <- (sb_cur$End - sb_cur$Start) / 60 - ((sb_cur$deepdive_n - 1)*5)
      plot_inputs[[count]] <- list(
        obs = iddi_dur[exp_ind], # minutes
        bas = iddi_dur[modified_baseline], # minutes
        info = info,
        what = "IDDI"
      )
      count <- count + 1
      
      # output prior IDDI
      if(exp_ind - 1 > 0) {
        plot_inputs[[count]] <- list(
          obs = iddi_dur[exp_ind - 1],
          bas = iddi_dur[modified_baseline],
          info = info,
          what = "prior_IDDI"
        )
        count <- count + 1
      }
      
      # output nbounce
      plot_inputs[[count]] <- list(
        obs = sb_cur$bounce_n[exp_ind], # n bounce dives in IDDI
        bas = sb_cur$bounce_n[modified_baseline], # n bounce dives in IDDI baseline
        info = info,
        what = "n_bounce"
      )
      count <- count + 1
      
      # output mean_bounce_duration
      plot_inputs[[count]] <- list(
        obs = sb_cur$average_bounce_nsamp[exp_ind] * 5, # minutes bounce dur
        bas = sb_cur$average_bounce_nsamp[modified_baseline] * 5, # minutes bounce dur baseline
        info = info,
        what = "mean_bounce_duration"
      )
      count <- count + 1
      
      # output mean_bounce_maxdepth
      plot_inputs[[count]] <- list(
        obs = sb_cur$bounce_meandepth[exp_ind], # meters
        bas = sb_cur$bounce_meandepth[modified_baseline], # meters baseline
        info = info,
        what = "mean_bounce_maxdepth"
      )
      count <- count + 1
      
      # output bounce_maxdepth
      plot_inputs[[count]] <- list(
        obs = sb_cur$bounce_maxdepth[exp_ind], # meters
        bas = sb_cur$bounce_maxdepth[modified_baseline], # meters baseline
        info = info,
        what = "bounce_maxdepth"
      )
      count <- count + 1
    }
  }
}

saveRDS(plot_inputs, file = "02_intermediate_outputs/plot_inputs.rds")
