set.seed(314)

library(dplyr)
library(readr)
library(stringr)

# load summary statistic fn.
source('sattag_summary.R')

# load additional components for simulation
source('segment_fn.R')

# Read the .rds file that define the template depth bins
template_bins <- readRDS(here::here("data/template_bins.rds"))
bin_starts = template_bins$center - template_bins$halfwidth

# load all window configurations
study_windows <- readRDS(here::here("data/study_windows.rds"))

# Load the tag information
# tag_info_all <- read_csv("data/tag_info_all-tags.csv", show_col_types = FALSE)
tag_indices_bline <- read_csv(here::here("data/tag_info_mult-segments_baseline.csv"), 
                              show_col_types = FALSE)
tag_indices_cee <- read_csv(here::here("data/sattag/tag_info_mult-segments_cee.csv"), 
                            show_col_types = FALSE)
tags_mult_bline_metadata <- read_csv(file = here::here('data/tag_info_mult-segments_baseline_metadata.csv'))

# Cull out only the files in tag_info
dive_dir <- "../../01_shared_data_products/filter_sattag/series/"
my_files <- data.frame(my_full_files = dir(dive_dir,
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE,
                                           full.names = TRUE),
                       my_file_names = dir(dive_dir,
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE))

my_files$my_files_deploy <- stringr::str_sub(my_files$my_file_names, start = 1, end = 8)  
idx <- which(my_files$my_files_deploy %in% tag_indices_bline$deployid)
my_files <- my_files[idx, ]

# Function to do depth labeling and enriching
depth_enrich = function(x, time_inds){
  
  d = lapply(time_inds, function(inds) {
    data.frame(depths = x[inds]) %>% 
      mutate(
        depth.bin = findInterval(x = depths, vec = bin_starts),
        ddepths = c(0, diff(depths)),
        ddepths.sign = sign(ddepths),
        descending = ddepths.sign > 0,
        ascending = ddepths.sign < 0,
        no_change = ddepths.sign == 0
      )  
  }
  )
  
  # segment dives
  # also a list - see above block
  # or put ALL of these steps into a giant lapply
  labels = lapply(d, function(d) {
    
    dive.segmentation(
      y = d$depth.bin, merge.ratio = .6, depth.bins = template_bins, 
      times = as.POSIXct(1:nrow(d), origin = '1970-01-01 00:00.00 UTC'), 
      timestep = 1)
    
  }
  )
  
  # compute dive types
  diveTypes = lapply(1:length(d), function(ind) {
    d[[ind]] %>% 
      mutate(diveId = labels[[ind]]) %>% 
      group_by(diveId) %>% 
      summarise(
        diveType = ifelse(max(depths) > 800, 'Deep', 'Shallow')
      )
  }
  )
  
  # merge dive types
  d = lapply(1:length(d), function(ind) {
    
    d[[ind]] %>% 
      mutate(diveId = labels[[ind]]) %>% 
      left_join(diveTypes[[ind]], by = 'diveId')
  }
  )
  
  return(d)
}

# function to analyze time series of depths
prepost = function(x, nwin, nlag, baseline_end, cee_idx, tag) {
  # Parameters:
  #  x - time series of depths
  #  nwin - number of timepoints in each pre/post window
  #  nlag - number of timepoints between each pre/post window
  #  baseline_end - list of indices of _last_ baseline observation(s)
  #  cee_idx - list of indices of each cee start; 
  # TODO: bring in the actual cee id value, rather than the integers
  # Possible TODO: but we might want to make cee_idx like
  # > tag_indices_cee
  # A tibble: 2 × 3
  # deployid potential_response_start_idx potential_response_end_idx
  # <chr>                 <dbl>            <dbl>
  #  1 CEE 01                  1              301
  #  2 CEE 02                589             1756
  # where potential_response_end_idx is the 24-hour rule, and if there's a gap in 
  # data collection, we cut this short(er): 589 = 301 + 24*12 assuming no gaps
  
  # Assemble list of baseline indices
  pre_inds <- vector(mode = 'list', length = nrow(cee_idx))
  for(j in 1:nrow(cee_idx)){
    
    pre_i <- seq(to = cee_idx$cee_start_idx[j] - 1, length.out = nwin)
    
    pre_inds[[j]] <- pre_i
  }
  
  
  
  # need to change this to loop over cee starts (which are TO BE read in separately), and right now it has
  # 1:1 correspondence with baseline indices. This well help with flexibility 
  # and gaps (which need to be processed before this function runs - hence elsewhere)
  # search for gap_after in the dsdive_fulltag repo for inspiration
  post_inds <- vector(mode = 'list', length = nrow(cee_idx))
  for(j in 1:nrow(cee_idx)){
    
    post_i <- seq(from = tail(unlist(cee_idx[j, 'cee_start_idx']), 1) + nlag + 1, length.out = nwin)
    
    # TODO: this is wrong
    # validate that enough data is present for analysis
    # if(max(c(pre_inds[[j]], post_i)) > length(x)) {
    #   stop('Not enough data to run analysis')
    # }
    
    # Assemble
    post_inds[[j]] <- post_i
  }
  
  # This will need to be a list for each baseline period
  # initial features
  # x will be subset to correspond to the start/end of this sequence of baseline obs
  # the baseline indices are portions of tag WITHOUT data gaps
  # have to assume that pre_inds has NO gaps
  
  ############################################################
  # Enrich Depth Data with Indices and Labels
  ############################################################
  d_pre <- depth_enrich(x, pre_inds)
  d_post <- depth_enrich(x, post_inds)
  
  # observed results
  # will need to be run for each of the cees; so a list to hold it
  observed_pre = lapply(d_pre, function(d) {
    sattag_summary(d)
  })
  
  observed_post = lapply(d_post, function(d) {
    sattag_summary(d)
  })
  
  # browser()
  
  # Assemble the null samples
  # pre_inds are the baseline periods
  null_samples_pre_list <- vector(mode = 'list', length = nrow(baseline_end))
  
  for(j in 1:nrow(baseline_end)){
    d <- depth_enrich(x, list(baseline_end$baseline_start_idx[j]:baseline_end$baseline_end_idx[j]))[[1]]
    
    null_samples_pre = do.call(rbind, lapply(1:nrow(d), 
                                             function(start) {
      # determine window inds
      pre_inds = seq(from = start, length.out = nwin)
      post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
      
      # skip processing if any indices are outside baseline window
      if(max(c(pre_inds,post_inds)) > nrow(d)) {
        return(NULL)
      }
      # compute summaries
      sattag_summary(d[pre_inds,])
    }))
    
    null_samples_pre_list[[j]] <- null_samples_pre
  }
  
  null_samples_pre <- do.call(rbind, null_samples_pre_list)
  
  
  # Post Windows
  # when it says null, we're only ever talking about baseline
  null_samples_post_list <- vector(mode = 'list', length = nrow(baseline_end))
  for(j in 1:nrow(baseline_end)){
    d <- depth_enrich(x, list(baseline_end$baseline_start_idx[j]:baseline_end$baseline_end_idx[j]))[[1]]
    
    null_samples_post = do.call(rbind, lapply(1:nrow(d), 
                                              function(start) {
    # determine window inds
    pre_inds = seq(from = start, length.out = nwin)
    post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
    # skip processing if any indices are outside baseline window
    if(max(c(pre_inds,post_inds)) > nrow(d)) {
      return(NULL)
    }
    # compute summaries
    sattag_summary(d[post_inds,])
  }))
    
    null_samples_post_list[[j]] <- null_samples_post
  }
  
  null_samples_post <- do.call(rbind, null_samples_post_list)
  
  ######################################################################
  # Make some plots
  # ind <- 3 # time deep diving
  # df <- data.frame(
  #   pre = null_samples_pre[,ind],
  #   post = null_samples_post[,ind]
  # )
  # 
  # pl = ggplot(df, aes(x = pre, y = post)) + 
  #   # overplotted raw data
  #   geom_point(col = 'grey', alpha = .3) +
  #   # bivariate density
  #   stat_density_2d(col = 'black') +
  #   # conditional baseline distribution used for test
  #   geom_vline(xintercept = observed_pre[[1]]$time_deep_diving, col = 'darkgreen') +# observed_pre[,ind]
  #   # pre/post pair observed during CEE
  #   geom_point(
  #     data = data.frame(
  #       pre = observed_pre[[1]]$time_deep_diving, #observed_pre[, ind], 
  #       post = observed_post[[1]]$time_deep_diving# observed_post[, ind]
  #     ),
  #     col = 'darkgreen'
  #   ) + 
  #   # axis labels and formatting
  #   xlab('Time Deep Diving - Pre')+ 
  #   ylab('Time Deep Diving - Post')+ 
  #   labs(subtitle = my_deployid)+
  #   theme_bw() + 
  #   theme(axis.title.y = element_text(angle = 0, vjust = .5)) + 
  #   coord_equal() 
  # 
  # f = file.path('output', 'scaled_source', 'kde_plots')
  # dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
  # ggsave(pl, 
  #        filename = file.path(f, paste('_bivariate_', tag, 
  #                                      '.pdf', sep ='')),
  #        width = 8, height = 8, dpi = 'print')
  ######################################################################
  
  #
  # compute tail probabilities
  #
  
  # This is going to be outside of a list; will process null_samples_pre
  # from all baseline periods together
  probs = do.call(rbind, lapply(1:ncol(null_samples_pre), function(stat_ind) {
    
    # ind indexes the statistic,
    # but we also need the cee_idx as well inside of this
    
    # bandwidths for kde, modified s.t. results are always non-zero if data
    # are not degenerate
    h = c(MASS::bandwidth.nrd(null_samples_pre[,stat_ind]),
          MASS::bandwidth.nrd(null_samples_post[,stat_ind]))
    if(h[1] == 0) {
      x = null_samples_pre[,stat_ind]
      h[1] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    if(h[2] == 0) {
      x = null_samples_post[,stat_ind]
      h[2] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    
    # bivariate kernel density estimate (KDE)
    dens = MASS::kde2d(
      x = null_samples_pre[,stat_ind], 
      y = null_samples_post[,stat_ind], 
      n = 500, 
      lims = c(range(c(null_samples_pre[,stat_ind], sapply(observed_pre, function(x){x[[stat_ind]]}) )), 
               range(c(null_samples_post[,stat_ind], sapply(observed_post, function(x){x[[stat_ind]]}) )),
      h = h
    ))
    
    # This line to end will be replicated for each CEE
    # wrap in a do.call(rbind, lapply(cee_ind) )
    # indexes to retrieve conditional distributions for KDE 
    kde_ind_pre = sapply(observed_pre, function(x) which.min(abs(dens$x - x[[stat_ind]])))
    kde_ind_post = sapply(observed_post, function(x) which.min(abs(dens$y - x[[stat_ind]]))) 
    
    # conditional distribution P(Post | Pre) and complement
    # Fix to account for multiple baselines with a single CEE
    if(nrow(cee_idx) > 1){
      C = rowSums(dens$z[kde_ind_pre,])
    } else {
      C = sum(dens$z[kde_ind_pre,])
    }
    
    cdf.conditional = t(sapply(kde_ind_pre, function(x) cumsum(dens$z[x,]) )) / C
    ccdf.conditional = t(sapply(kde_ind_pre, function(x) rev(cumsum(rev(dens$z[x,]))) )) / C
    
    do.call(rbind, lapply(1:length(kde_ind_post), 
                          function(cee_idx) 
                          {
                            data.frame(
                            p.left = cdf.conditional[cee_idx, kde_ind_post[cee_idx]],
                            p.right = ccdf.conditional[cee_idx, kde_ind_post[cee_idx]],
                            stat = factor(colnames(null_samples_pre)[stat_ind]),
                            cee_idx = cee_idx
                            )
                          }  
                          ))
    
  }))
  
  probs
}

f = file.path('output', 'scaled_source_24hr')
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

# Loop over animals - Proper function looping; 
# The for loop below is more making the plot for the scaled source ms
for(i in 1:nrow(my_files)){

# Testing for plots
# my_files <- my_files[which(my_files$my_files_deploy %in% c("ZcTag093", "ZcTag096")), ]
# for(i in 1:nrow(my_files)){

  animal <- read_csv(my_files$my_full_files[i], show_col_types = FALSE)
  my_deployid <- str_extract(unique(animal$DeployID), 'ZcTag[0-9]+')
  if(my_deployid == "ZcTag102") next()
  if(my_deployid == "ZcTag073"){
    
    zctag073_Day <- "18-Aug-2018"
    zctag073_Time <- "23:55:00"
    idx <- which(animal$Day == zctag073_Day & as.character(animal$Time) == zctag073_Time)
    animal <- animal[1:idx, ]
  }
  depths <- animal$Depth
  
  baseline_end <- tag_indices_bline %>% 
    filter(deployid == my_deployid)
  cee_idx <- tag_indices_cee %>% 
    filter(deployid == my_deployid)
  
  # Test for gaps and split up baseline end
  gaps <- which(diff(animal$Date) > 300)
  
  # check and remove any extra-baseline gaps
  gapscheck <- list()
  for(q in 1:nrow(baseline_end)) {
    gapscheck[[q]] <- dplyr::between(gaps, baseline_end$baseline_start_idx[q], baseline_end$baseline_end_idx[q])
  }
  
  gapscheck <- do.call('rbind', gapscheck)
  gaps <- gaps[apply(gapscheck, 2, any)]
  
  # work over gaps if they occur during the baseline
  if(length(gaps)>0){
    
    outlist <- vector(mode = "list", length = length(gaps))
    new_baseline_df <- baseline_end
    for (j in seq_along(gaps)) {
      column_name <- paste0("gap_", j)
      baseline_end_df <- new_baseline_df %>%
        mutate(!!column_name := gaps[j] >= baseline_start_idx & gaps[j] <= baseline_end_idx) 

      # row to change
      int_vec_df <- baseline_end_df %>% 
        filter(if_all(starts_with("gap")))
      
      # row to keep
      df_keep <- baseline_end_df %>% 
        filter(!if_all(starts_with("gap"))) %>% 
        select(-starts_with("gap"))
        
      # Split out the vector based on the gaps
      split_vectors <- vector(mode = 'list', length = 2)
      start_index <- int_vec_df$baseline_start_idx
      end_index <- int_vec_df$baseline_end_idx
      gap_end_index <- gaps[j]
      split_vectors[[1]] <- c(start_index, gap_end_index)
      split_vectors[[2]] <- c((gap_end_index + 1), end_index)
      
      # Assemble the interim updated row(s)
      int_vec_out <- data.frame(deployid = my_deployid,
                                baseline_start_idx = c(split_vectors[[1]][1], split_vectors[[2]][1]),
                                baseline_end_idx = c(split_vectors[[1]][2], split_vectors[[2]][2]))
      
      # Assemble the updated data frame
      new_baseline_df <- bind_rows(df_keep, int_vec_out) %>% 
        arrange(baseline_start_idx)
      outlist[[j]] <- new_baseline_df
      
      } # end seq along gaps
    
    baseline_end <- new_baseline_df
    
    } # end working over the gaps 


  
  # Run pre-post function
  tag_results <-  do.call(rbind, lapply(1:nrow(study_windows), function(wid) {
    
    # print(wid)
    # extract window config
    wlen = study_windows$window_length[wid]
    wlag = study_windows$response_lag[wid]
    # convert to number of observations
    nwin = wlen * 60 / 5
    nlag = wlag * 60 / 5
    # analyze datasets
    cbind(
      # simulation details
      tag = my_deployid,
      window_len = wlen,
      lag_len = wlag,
      prepost(x = depths, nwin = nwin, nlag = nlag, baseline_end, cee_idx, tag = my_deployid) #baseline_end, cee_idx
    )
  }))
  
  fname = file.path(f, paste('samples_', my_deployid, '.rds', sep = ''))
  saveRDS(tag_results, file = fname)
  
}
