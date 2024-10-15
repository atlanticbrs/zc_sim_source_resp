###
# setup_definedives.r
# load metadata and re-factor dive data (beh + ser) 



### constants
SERIES_MAXIMUM_DEPTH_DIVE_ST_OR_EN <- 200 # arbitrary, but works well for this dataset

DEEP_DIVE_DEFINITION_METERS  <- 800 # based on shearer et al. 2019 which is based on baird's observations
LONG_DIVE_DEFINITION_MINUTES <-  33 # based on shearer et al. 2019 which is based on baird's observations
LONG_DIVE_DEFINITION_SAMPLES <-   6 # based on 5 minute samples should be 6
BEHAVIOR_MESSAGE_GAP_SECONDS <- 120 # 2 minutes because each message can be as much as 1 min off


### import data and make dates
#cee metadata import
cee <- read.table("../../00_data_input/cee/cee_metadata_flat.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
scaled <- cee[cee$cee_type == "simulated mfas" | cee$cee_type == "control", ]

# add datenum
scaled$st_datenum <- as.numeric(as.POSIXct(paste(scaled$date_YYYYMMDD, scaled$start_HHMMSS), format = "%Y%m%d %H%M%S", tz = "UTC"))
scaled$en_datenum <- as.numeric(as.POSIXct(paste(scaled$date_YYYYMMDD, scaled$end_HHMMSS), format = "%Y%m%d %H%M%S", tz = "UTC"))

# behavior data import
dpath <- "../../01_shared_data_products/filter_sattag/behavior/"
fnames <- list.files(dpath)
fnames <- fnames[grepl("-Behavior\\.", fnames)] # just the behavior

beh <- list()
for(i in 1:length(fnames)) {
  beh[[i]] <- read.table(file.path(dpath, fnames[i]), header = TRUE, sep = ',', stringsAsFactors = FALSE)
}

# series data import
dpath <- "../../01_shared_data_products/filter_sattag/series/"
fnames <- list.files(dpath)
fnames <- fnames[grepl("-Series\\.", fnames)] # just the series

ser <- list()
for(i in 1:length(fnames)) {
  ser[[i]] <- read.table(file.path(dpath, fnames[i]), header = TRUE, sep = ',', stringsAsFactors = FALSE)
}



### series dive definitions
# simple algorithm looking for peaks at 200 m and above
series_dive_definitions <- list()
pb <- txtProgressBar(style = 3)
for(q in 1:length(ser)) {
setTxtProgressBar(pb, q/length(ser))
  s1 <- ser[[q]]
  
  # define data gaps
  s1$tdiff <- c(0, s1$Date[2:nrow(s1)] - s1$Date[1:(nrow(s1) - 1)])
  s1$stretchid <- cumsum(s1$tdiff > 300) + 1
  
  # sign of the first derivative gives inflection points
  dts1 <- c(NA, diff(-s1$Depth))
  dts1_sign <- sign(dts1)
  s1$sign <- dts1_sign
  s1$dd <- abs(dts1)
  
  working_range <- 2:(nrow(s1) - 1)
  
  # any of these are acceptable transitions
  # that might be the start or end of a dive
  pass1 <- dts1_sign[working_range] == 1 & dts1_sign[working_range + 1] == -1
  pass2 <- dts1_sign[working_range] == 1 & dts1_sign[working_range + 1] == 0
  pass3 <- dts1_sign[working_range] == 0 & dts1_sign[working_range + 1] == -1
  
  # check the conditions for transitions and not over a gap
  pass_all <- pass1 | pass2 | pass3
  
  # most likely a change in the dive if they are also shallow
  pass_final <- pass_all & s1$Depth[working_range] <= SERIES_MAXIMUM_DEPTH_DIVE_ST_OR_EN
  pass_final <- c(FALSE, pass_final, FALSE)
  
  # define an index based on these dives
  diveid <- cumsum(pass_final)
  udiveid <- unique(diveid)
  ndiveid <- length(udiveid)
  
  # add back in the dataframe
  s1$pass_final <- pass_final
  s1$diveid <- diveid
  
  # drop each dive into an element in a list
  # skip the first one and the last one because could be in the middle
  count <- 1
  dives <- list()
  curstretchid <- 1
  for(i in 2:(ndiveid-1)) {
    dese <- which(diveid == udiveid[[i]])
    dese <- c(dese, dese[length(dese)] + 1)
    
    # check to make sure we haven't started a new stretch (crossed a gap)
    # if yes then strike this dive because it might be in the middle of a dive
    if(max(s1$stretchid[dese]) == curstretchid) {
      dives[[count]] <- s1[dese, ]
      dives[[count]]$diveid <- count
      count <- count + 1
    } else {
      curstretchid <- max(s1$stretchid[dese])
    }
  }
  
  # output
  series_dive_definitions[[q]] <- list(DeployID = s1$DeployID[1], dives = dives)
}
close(pb)



### calculate some dive by dive metrics for series
series_dive_summary <- list()

pb <- txtProgressBar(style = 3)
for(i in 1:length(series_dive_definitions)) {
setTxtProgressBar(pb, i/length(series_dive_definitions))
  curid <- series_dive_definitions[[i]]$DeployID
  curdives <- series_dive_definitions[[i]]$dives
  
  dat_tmp <- data.frame(DeployID = rep(curid, length(curdives)))
  dat_tmp$Start <- sapply(curdives, function(x) min(x$Date))
  dat_tmp$End <- sapply(curdives, function(x) max(x$Date))
  dat_tmp$nsamp <- sapply(curdives, nrow)
  dat_tmp$maxdepth <- sapply(curdives, function(x) max(x$Depth))
  dat_tmp$nascent <- sapply(curdives, function(x) {
    md <- max(which(max(x$Depth) == x$Depth))
    length(md:nrow(x))
  })
  dat_tmp$ndecent <- sapply(curdives, function(x) {
    md <- which(max(x$Depth) == x$Depth)[1]
    length(1:md)
  })
  
  # find gaps
  dat_tmp$tdiff <- c(0, dat_tmp$Start[2:nrow(dat_tmp)] - dat_tmp$End[1:(nrow(dat_tmp) - 1)])
  dat_tmp$stretchid <- cumsum(dat_tmp$tdiff != 0) + 1
  
  # output
  series_dive_summary[[i]] <- dat_tmp
}
close(pb)



### summarize characteristics of bouts for series
series_bout_summary <- list()
pb <- txtProgressBar(style = 3)
for(i in 1:length(series_dive_summary)) {
setTxtProgressBar(pb, i/length(series_dive_summary))
  curdivesum <- series_dive_summary[[i]]
  
  # only want to define 'deep' dives as deep and long
  curdivesum$deep <- (curdivesum$maxdepth > DEEP_DIVE_DEFINITION_METERS) & (curdivesum$nsamp > LONG_DIVE_DEFINITION_SAMPLES)
  curdivesum$boutid <- cumsum(curdivesum$deep)
  
  # truncate the first and last because we don't know
  curdivesum <- curdivesum[curdivesum$boutid > 0 & curdivesum$boutid < max(curdivesum$boutid), ]
  
  # make a dataframe for the bouts
  ubouts <- unique(curdivesum$boutid)
  nbouts <- length(ubouts)
  bout_summary_tmp <- data.frame(DeployID = rep(curdivesum$DeployID[1], nbouts), boutid = ubouts)
  
  deepdives <- curdivesum[curdivesum$deep, ]
  bouts_tmp <- split(curdivesum, curdivesum$boutid)
  
  
  bout_summary_tmp$hasgap <- sapply(bouts_tmp, function(x) length(unique(x$stretchid)) > 1)
  bout_summary_tmp$Start <- sapply(bouts_tmp, function(x) min(x$Start))
  bout_summary_tmp$End <- sapply(bouts_tmp, function(x) max(x$End))
  
  bout_summary_tmp$deepdive_depth <- deepdives$maxdepth
  bout_summary_tmp$deepdive_n <- deepdives$nsamp
  bout_summary_tmp$bounce_n <- sapply(bouts_tmp, function(x) {
    length(which(!x$deep & x$maxdepth > 150))
  })
  
  bout_summary_tmp$average_bounce_nsamp <- sapply(bouts_tmp, function(x) {
    mean(x$nsamp[!x$deep & x$maxdepth > 150])
  })
  bout_summary_tmp$bounce_maxdepth <- sapply(bouts_tmp, function(x) {
    tmpout <- max(x$maxdepth[!x$deep])
    if(is.infinite(tmpout)) tmpout <- NA
    tmpout
  })
  bout_summary_tmp$bounce_meandepth <- sapply(bouts_tmp, function(x) {
    mean(x$maxdepth[!x$deep & x$maxdepth > 150])
  })
  
  # remove any bouts that have a gap in them
  bout_summary_tmp <- bout_summary_tmp[!bout_summary_tmp$hasgap, ]
  
  series_bout_summary[[i]] <- bout_summary_tmp
}
close(pb)



### refactor behavior into a deep dive > DEEP_DIVE_DEFINITION_METERS (> 800 m) and IDDI format
behavior_bout_summary <- list()
pb <- txtProgressBar(style = 3)
for(i in 1:length(beh)) {
setTxtProgressBar(pb, i/length(beh))
  b1 <- beh[[i]]
  b1 <- b1[b1$What != "Message", ]
  
  # average depth and duration
  b1$Depth <- (b1$DepthMin + b1$DepthMax) / 2
  b1$Duration <- b1$End - b1$Start
  
  ### find the gaps
  # going to use a threshold of 2 minutes just for some wiggle room as
  # the clock is imprecise between message blocks
  b1$tdiff <- c(0, b1$Start[2:nrow(b1)] - b1$End[1:(nrow(b1) - 1)])
  b1$stretchid <- cumsum(b1$tdiff > BEHAVIOR_MESSAGE_GAP_SECONDS) + 1
  
  # split by stretchid (gaps)
  bl <- split(b1, b1$stretchid)
  
  bouts_list <- list()
  for(q in 1:length(bl)) {
    b2 <- bl[[q]]
  
    # look for any dives under DEEP_DIVE_DEFINITION_METERS (800) that slipped through
    # and convert them to surface
    b2$What[b2$Depth <= DEEP_DIVE_DEFINITION_METERS] <- "Surface"
    
    # look for any dives under LONG_DIVE_DEFINITION_MINUTES (33) that slipped through
    # and convert them to surface
    b2$What[b2$Duration <= LONG_DIVE_DEFINITION_MINUTES] <- "Surface"
    
    # look for adjacent surfaces to combine
    surfid <- cumsum(b2$What == "Dive")
    surfid[b2$What == "Dive"] <- NA
    surfid_table <- table(surfid)
    usurfid <- as.numeric(names(surfid_table))[surfid_table > 1]
  
    if(any(surfid > 0)) {
      b2$keep <- TRUE
  
      for(cursurfid in usurfid) {
        dese <- which(surfid == cursurfid)
      
        # calculate new end and start
        st <- min(b2$Start[dese])
        en <- max(b2$End[dese])
      
        # assign new end and start and duration
        b2$Start[dese[1]] <- st
        b2$End[dese[1]] <- en
        b2$Duration[dese[1]] <- en - st
      
        # designate obsolete for removal
        b2$keep[dese[2:length(dese)]] <- FALSE
      }
    
      # remove obsolete surfaces
      b2 <- b2[b2$keep, ]
      
      # if the first entry is a surface kill it because we might be in the middle
      if(b2$What[1] == 'Surface') b2 <- b2[-1, ]
    
      ### make a new data frame where each row is a dive and the following IDDI
      dive <- b2[b2$What == "Dive", ]
      surf <- b2[b2$What == "Surface", ]
  
      bouts <- data.frame(
        DeployID = dive$DeployID,
        dive_st = dive$Start,
        dive_en = dive$End,
        shape = dive$Shape,
        depth = dive$Depth,
        dive_duration = dive$Duration,
        stretchid = dive$stretchid,
        surf_st = NA,
        surf_en = NA,
        surf_duration = NA
      )
      
      bouts$surf_st[1:nrow(surf)] <- surf$Start
      bouts$surf_en[1:nrow(surf)] <- surf$End
      bouts$surf_duration <- bouts$surf_en - bouts$surf_st
      
      bouts_list[[q]] <- bouts
    } # is this close brace in the right place? is surfid ever 0?
  }
    
  # save output
  bouts_df <- do.call('rbind', bouts_list)
  behavior_bout_summary[[i]] <- bouts_df
}
close(pb)



### output intermediate files
save(
  scaled, ser, beh, # input data
  series_dive_definitions, # intermediate file for error checking
  series_dive_summary, series_bout_summary, # series output for calculation
  behavior_bout_summary, # behavioroutput for calculation
  
  file = "02_intermediate_outputs/setup_define_out.rdata" # output file
)
