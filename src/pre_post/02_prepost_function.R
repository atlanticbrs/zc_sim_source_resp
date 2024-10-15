set.seed(314)

library(tidyverse)
library(dplyr)

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
tag_info <- read_csv("data/tag_info.csv", show_col_types = FALSE)

# Cull out only the files in tag_info
my_files <- data.frame(my_full_files = dir(here::here("../../01_shared_data_products/filter_sattag/series/"),
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE,
                                           full.names = TRUE),
                       my_file_names = dir(here::here("../../01_shared_data_products/filter_sattag/series/"),
                                           pattern = glob2rx("Zc*Series.csv"),
                                           recursive = TRUE))

my_files$my_files_deploy <- stringr::str_sub(my_files$my_file_names, start = 1, end = 8)  
idx <- which(my_files$my_files_deploy %in% tag_info$deployid)
my_files <- my_files[idx, ]



# function to analyze time series of depths
prepost = function(x, nwin, nlag, baseline_end) {
  # Parameters:
  #  x - time series of depths
  #  nwin - number of timepoints in each pre/post window
  #  nlag - number of timepoints between each pre/post window
  #  baseline_end - index of last baseline observation
  
  # TODO: make this a list for each segment
  # observed pre/post
  pre_inds = seq(to = baseline_end, length.out = nwin)
  post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
  
  # validate that enough data is present for analysis
  if(max(c(pre_inds,post_inds)) > length(x)) {
    stop('Not enough data to run analysis')
  }
  
  #
  # enrich time series
  #
  
  # This will need to be a list for each baseline period
  # initial features
  # x will be subset to correspond to the start/end of this sequence of baseline obs
  d = data.frame(depths = x) %>% 
    mutate(
      depth.bin = findInterval(x = depths, vec = bin_starts),
      ddepths = c(0, diff(depths)),
      ddepths.sign = sign(ddepths),
      descending = ddepths.sign > 0,
      ascending = ddepths.sign < 0,
      no_change = ddepths.sign == 0
    )
  
  # segment dives
  # also a list
  labels = dive.segmentation(
    y = d$depth.bin, merge.ratio = .6, depth.bins = template_bins, 
    times = as.POSIXct(1:nrow(d), origin = '1970-01-01 00:00.00 UTC'), 
    timestep = 1
  )
  
  # compute dive types
  diveTypes = d %>% 
    mutate(diveId = labels) %>% 
    group_by(diveId) %>% 
    summarise(
      diveType = ifelse(max(depths) > 800, 'Deep', 'Shallow')
    )
  
  # merge dive types
  d = d %>% 
    mutate(diveId = labels) %>% 
    left_join(diveTypes, by = 'diveId')
  
  #
  # compute window summaries
  #
  
  # observed results
  # will need to be run for each of the cees; so a list to hold it
  observed_pre = sattag_summary(d[pre_inds, ])
  observed_post = sattag_summary(d[post_inds, ])
  
  # This will get run for every baseline segment
  # but won't be stored as a list; will be a single df
  # null distribution samples
  null_samples_pre = do.call(rbind, lapply(1:baseline_end, function(start) {
    # determine window inds
    pre_inds = seq(from = start, length.out = nwin)
    post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
    # skip processing if any indices are outside baseline window
    if(max(c(pre_inds,post_inds)) > bend) {
      return(NULL)
    }
    # compute summaries
    sattag_summary(d[pre_inds,])
  }))
  null_samples_post = do.call(rbind, lapply(1:baseline_end, function(start) {
    # determine window inds
    pre_inds = seq(from = start, length.out = nwin)
    post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
    # skip processing if any indices are outside baseline window
    if(max(c(pre_inds,post_inds)) > bend) {
      return(NULL)
    }
    # compute summaries
    sattag_summary(d[post_inds,])
  }))
  
  ######################################################################
  # Make some plots
  ind <- 3 # time deep diving
  df <- data.frame(
    pre = null_samples_pre[,ind],
    post = null_samples_post[,ind]
  )

  pl = ggplot(df, aes(x = pre, y = post)) +
    # overplotted raw data
    geom_point(col = 'grey', alpha = .3) +
    # bivariate density
    stat_density_2d(col = 'black') +
    # conditional baseline distribution used for test
    geom_vline(xintercept = observed_pre[,ind], col = 'darkgreen') +# observed_pre[,ind]
    # pre/post pair observed during CEE
    geom_point(
      data = data.frame(
        pre = observed_pre[, ind],
        post = observed_post[, ind]
      ),
      col = 'darkgreen'
    ) +
    # axis labels and formatting
    xlab('Time Deep Diving - Pre')+
    ylab('Time Deep Diving - Post')+
    labs(subtitle = tag)+
    theme_bw() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5)) +
    coord_equal()

  f = file.path('output', 'scaled_source', 'kde_plots')
  dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
  ggsave(pl,
         filename = file.path(f, paste('_bivariate_', tag,
                                       '.pdf', sep ='')),
         width = 8, height = 8, dpi = 'print')
  ######################################################################
  
  #
  # compute tail probabilities
  #
  
  probs = do.call(rbind, lapply(1:ncol(null_samples_pre), function(ind) {
    
    # bandwidths for kde, modified s.t. results are always non-zero if data
    # are not degenerate
    h = c(MASS::bandwidth.nrd(null_samples_pre[,ind]),
          MASS::bandwidth.nrd(null_samples_post[,ind]))
    if(h[1] == 0) {
      x = null_samples_pre[,ind]
      h[1] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    if(h[2] == 0) {
      x = null_samples_post[,ind]
      h[2] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    
    # bivariate kernel density estimate (KDE)
    dens = MASS::kde2d(
      x = null_samples_pre[,ind], 
      y = null_samples_post[,ind], 
      n = 500, 
      lims = c(range(c(null_samples_pre[,ind], observed_pre[[ind]])), 
               range(c(null_samples_post[,ind], observed_post[[ind]]))),
      h = h
    )
    
    # This line to end will be replicated for each CEE
    # wrap in a do.call(rbind, lapply(cee) )
    # indexes to retrieve conditional distributions for KDE 
    kde_ind_pre = which.min(abs(dens$x - observed_pre[[ind]]))
    kde_ind_post = which.min(abs(dens$y - observed_post[[ind]]))
    
    # conditional distribution P(Post | Pre) and complement
    C = sum(dens$z[kde_ind_pre,])
    cdf.conditional = cumsum(dens$z[kde_ind_pre,]) / C
    ccdf.conditional = rev(cumsum(rev(dens$z[kde_ind_pre,]))) / C
    
    # This will need to be annotated with cee name
    # bivariate p-values
    data.frame(
      p.left = cdf.conditional[kde_ind_post],
      p.right = ccdf.conditional[kde_ind_post],
      stat = factor(colnames(null_samples_pre)[ind])
    )
  }))
  
  probs
}


f = file.path('output', 'scaled_source')
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

# Loop over animals
# for(i in 1:nrow(my_files)){
  
my_files <- my_files[which(my_files$my_files_deploy %in% c("ZcTag093", "ZcTag096")), ]
for(i in 1:nrow(my_files)){
  
  tag = my_files$my_files_deploy[i]
  stat = 'time_deep_diving'
  
  animal <- read_csv(my_files$my_full_files[i], show_col_types = FALSE)
  depths <- animal$Depth
  
  bend <- tag_info %>%
    filter(deployid == gsub("_DUML", "", unique(animal$DeployID))) %>%
    pull(cee_start_idx)

  deployid <- tag_info %>%
    filter(deployid == gsub("_DUML", "", unique(animal$DeployID))) %>%
    pull(deployid)
  
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
      tag = deployid,
      window_len = wlen,
      lag_len = wlag,
      prepost(x = depths, nwin = nwin, nlag = nlag, baseline_end = bend)
    )
  }))
  
  fname = file.path(f, paste('samples_', deployid, '.rds', sep = ''))
  saveRDS(tag_results, file = fname)
  
}
