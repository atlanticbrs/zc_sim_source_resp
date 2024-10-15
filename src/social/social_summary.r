###
# social_summary.r
# take a look at group dynamics around cees



### helper functions
# plot pretty (a), (b), ..., etc. legends on par multipanel plots
addsubplotleg <- function(lab, cex = 1.5, ...) {
  devin <- dev.size("in")
  x <- grconvertX(c(0, devin[1]), from = "in", to = "user")
  y <- grconvertY(c(0, devin[2]), from = "in", to = "user")
  
  fig <- par('fig')
  x <- x[1] + (x[2] - x[1]) * fig[1:2]
  y <- y[1] + (y[2] - y[1]) * fig[3:4]
  
  x <- x[1] + strwidth(lab, cex = cex)
  y <- y[2] - strheight(lab, cex = cex)
  text(x, y, lab, cex = cex, xpd = NA, ...)
}



### load data
cee <- read.table("../../00_data_input/cee/cee_metadata_flat.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
tagmeta <- read.table("../../01_shared_data_products/filter_sattag/metadata.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
rls <- read.table("../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

key <- read.table("00_data/grp_sightings_before_during_after.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
sight <- read.table("00_data/dumldb_indivsighting_prep.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
agesex <- read.table("00_data/zca_sex_20240103_DMW.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# make an rl key
rlkey <- paste(rls$deployid, rls$cee_id)
rlkey <- sub("_DUML", "", sub("Zc", "", rlkey))
 
# make sighting key and format date
sight$date_posix <- as.POSIXct(sight$Date, tz = "UTC", format = "%m/%d/%Y")
sight$key <- paste(format(sight$date_posix, "%Y%m%d"), toupper(sight$Sightno), sep = "_")

# formalize id names
key$Photo.ID <- paste0(key$Photo.ID, "_HAT")

# save matrices
info <- list()

for(i in 1:nrow(key)) {
  cur <- key[i, ]
  
  sightlist <- strsplit(c(cur$PRE, cur$CEE, cur$AFTER), " ")
  sightvec <- unlist(sightlist)
  
  ids <- list()
  for(s in 1:length(sightvec)) {
    ids[[s]] <- sight$id_code[sight$key == sightvec[s]]
  }
  
  allids <- sort(unique(unlist(ids)))
  
  sightmat <- matrix(0, length(allids), length(sightvec))
  colnames(sightmat) <- sightvec
  rownames(sightmat) <- allids
  
  # to highlight row and col of focal and the exposure
  exposuremat <- sightmat
  exposuremat[, which(colnames(exposuremat) %in% sightlist[[2]])] <- 1
  
  focalmat <- sightmat
  focalmat[which(rownames(focalmat) == key$Photo.ID[i]), ] <- 1

  # populate sightmat
  for(s in 1:length(ids)) {
    curids <- ids[[s]]
    
    sightmat[curids, s] <- 1
  }
  
  # only take the max count sighting from each day
  # if sighting names are the same it takes the first
  # this is not sorted by time though
  sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
  
  dese <- unlist(lapply(lapply(split(colSums(sightmat), sight_days), which.max), names))
  sightmat <- sightmat[, dese]
  focalmat <- focalmat[, dese]
  
  # set to NA anything that has no sightings in it (only would be a cee without obs)
  sightmat[, colSums(sightmat) == 0] <- NA
  
  # calculate lags
  # only include sightings that are within +/- 14 days
  sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
  dates <- as.Date(sight_days, format = "%Y%m%d")
  cee_date <- as.Date(as.character(cur$CEE.Date.1), format = "%Y%m%d")
  lags <- dates - cee_date
  dese <- abs(lags) <= 14
  
  sightmat <- sightmat[, dese]
  focalmat <- focalmat[, dese]
  
  # if a row doesn't have any sightings remove it
  focalmat <- focalmat[rowSums(sightmat, na.rm = TRUE) != 0, ]
  sightmat <- sightmat[rowSums(sightmat, na.rm = TRUE) != 0, ]
  
  xx <- seq(0, 1, length = ncol(sightmat))
  yy <- seq(0, 1, length = nrow(sightmat))
  
  # make the rownames labels shorter
  rnames <- rownames(sightmat)
  rnames <- sub("Zca_", "", rnames)
  rnames <- sub("_HAT", "", rnames)
  rnames[grep("Rakes", rnames)] <- paste0(sub("Rakes_", "", rnames[grep("Rakes", rnames)]), 'R')
  rnames <- sapply(strsplit(rnames, "_"), function(x) x[length(x)])
  rnames <- sub("Number", "no", rnames)
  
  # match any tag ids
  rownames_labs <- rnames
  rownames_labs <- sub("NA ", "", rownames_labs)
  
  # match agesex info
  matches <- match(rownames(sightmat), paste0(agesex$ID, "_HAT"))
  agesex_lab <- paste0(agesex$Age.Class[matches], agesex$Sex[matches])
  agesex_lab <- gsub("NA", "", agesex_lab)
  agesex_lab <- gsub(" ", "", agesex_lab)
  agesex_lab <- gsub("Unknown", "U", agesex_lab)
  agesex_lab <- sub("Adult", "A", agesex_lab)
  agesex_lab <- sub("Male", "M", agesex_lab)
  agesex_lab <- sub("Female", "F", agesex_lab)
  
  rownames_labs <- paste(rownames_labs, agesex_lab)
  
  # calculate +/- days from cee
  sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
  dates <- as.Date(sight_days, format = "%Y%m%d")
  cee_date <- as.Date(as.character(cur$CEE.Date.1), format = "%Y%m%d")
  lags <- dates - cee_date
  
  # calculate percent change
  # add a faux col to the beginning where no one is sighting
  sightmat <- sightmat[, colSums(apply(sightmat, 2, is.na)) == 0]
  sightmat_faux <- cbind(-1, sightmat)
  
  if(ncol(sightmat_faux) > 2) {
    sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
    dates <- as.Date(sight_days, format = "%Y%m%d")
    cee_date <- as.Date(as.character(cur$CEE.Date.1), format = "%Y%m%d")
    lags <- dates - cee_date
      
    diffmat <- sightmat_faux[, 2:ncol(sightmat_faux)] != sightmat_faux[, 1:(ncol(sightmat_faux) - 1 )]
    grpsize <- colSums(sightmat)
    changes <- colSums(diffmat)
    time_la <- as.numeric(lags - min(lags))
    totchan <- sum(changes[2:length(changes)]) 
    pchange <- changes / time_la / totchan / (ncol(diffmat) - 1)
    pchange[1] <- NA
    names(pchange) <- lags
    
    changes[1] <- NA
    names(changes) <- lags
    
    info[[i]] <- list(
      pchange = pchange,
      rawchange = changes,
      grpsize = grpsize,
      cee_date = as.character(format(cee_date, "%Y%m%d")),
      id = cur$Tag.ID,
      sightmat = sightmat,
      focalmat = focalmat,
      sightlist = sightlist,
      i = i,
      xx = xx,
      yy = yy,
      cur = cur,
      rownames_labs = rownames_labs
    )
  }
}

# prepare for plotting
info <- info[!sapply(info, is.null)]
sightmats <- lapply(info, '[[', 'rawchange')
grpsizes  <- lapply(info, '[[', 'grpsize')


lags <- lapply(sightmats, function(x) as.numeric(names(x)))
lags.old <- lags

# this is we want to squish the days without reference to
# actual time elapsed
lags <- lapply(lags, function(x) {
  before <- x[x < 0]
  before <- rev(seq(-1, -length(before), by = -1))
  after <- x[x > 0]
  if(length(after) > 0) {
    after <- seq(1, length(after), by = 1)
  } else {
    after <- NULL
  }
  
  zero <- NULL
  if(length(which(x == 0)) != 0) zero <- 0
  
  c(before, zero, after)  
})

# lets use actual time for now
lags <- lags.old

alllags <- seq(min(unlist(lags)), max(unlist(lags)), by = 1)

bigmat <- matrix(NA, nrow = length(sightmats), ncol = length(alllags))
grpmat <- bigmat

colnames(bigmat) <- alllags
colnames(grpmat) <- alllags

for(i in 1:length(sightmats)) {
  cur <- sightmats[[i]]
  grp <- grpsizes[[i]]
  
  bigmat[i, as.character(lags[[i]])] <- cur
  grpmat[i, as.character(lags[[i]])] <- grp
}
 
bigmat <- bigmat[-3, ]
grpmat <- grpmat[-3, ]
info <- info[-3]


# add in cee labels
ceeid <- cee$cee_id[match(sapply(info, '[[', 'cee_date'), cee$date_YYYYMMDD)]
ceetype <- cee$cee_type[match(sapply(info, '[[', 'cee_date'), cee$date_YYYYMMDD)]
deployid  <- sapply(info, '[[', 'id')
deployid <- sub("Zc", "", deployid)
deployid <- sub("^_", "", deployid)

matches <- match(paste(deployid, ceeid), rlkey)
rl_bins <- rls$rl_bin[matches]
rl_bins <- factor(rl_bins, levels = c("control", "<100", "100-120", ">120"))

rownames(bigmat) <- paste(deployid, ceeid, paste0("(", rl_bins, ")"))

oo <- order(rl_bins)
bigmat <- bigmat[oo, ]
grpmat <- grpmat[oo, ]



### social summary plot
pdf("social_summary.pdf", width = 10, height = 8)

lay <- layout(matrix(c(rep(1, 4), rep(2, 2), rep(3, 2)), nrow = 2, ncol = 4, byrow = TRUE))
par(mar = c(4.1, 15.1, 1.1, 1.1), family = "serif")



### big matrix
image(t(bigmat), axes = FALSE, col = hcl.colors(10, "blues", rev = TRUE))

# add in grey boxes for the NAs
greymat <- grpmat*0
greymat[is.na(grpmat)] <- 1
image(t(greymat), add = TRUE, col = c(NA, "grey75"))


# fix names
rownames(bigmat) <- gsub("Tag", "STag", rownames(bigmat))
rownames(bigmat) <- sub("19_218a", "DTag19_218", rownames(bigmat))
rownames(bigmat) <- sub("20_232a", "DTag20_232", rownames(bigmat))


rownames(bigmat)[rev(grep(">120", rownames(bigmat)))[1]] <- paste(      ">120 dB SPL        ", rownames(bigmat)[rev(grep(">120", rownames(bigmat)))[1]])
rownames(bigmat)[rev(grep("100-120", rownames(bigmat)))[1]] <- paste("100-120 dB SPL        ", rownames(bigmat)[rev(grep("100-120", rownames(bigmat)))[1]])
rownames(bigmat)[rev(grep("<100", rownames(bigmat)))[1]] <- paste(      "<100 dB SPL        ", rownames(bigmat)[rev(grep("<100", rownames(bigmat)))[1]])
rownames(bigmat)[rev(grep("control", rownames(bigmat)))[1]] <- paste(       "control        ", rownames(bigmat)[rev(grep("control", rownames(bigmat)))[1]])

rownames(bigmat) <- gsub(".\\(>120\\)$", "", rownames(bigmat))
rownames(bigmat) <- gsub(".\\(100-120\\)$", "", rownames(bigmat))
rownames(bigmat) <- gsub(".\\(control\\)$", "", rownames(bigmat))

modifiers <- rep("", length(colnames(bigmat)))
modifiers[sign(as.numeric(colnames(bigmat))) == 1] <- "+"
laglabels <- paste0(modifiers, colnames(bigmat))

box()
axis(1, at = seq(0, 1, length = ncol(bigmat)), lab = laglabels)
axis(2, at = seq(0, 1, length = nrow(bigmat)), lab = rownames(bigmat), las  = 1)
usr <- par('usr')
abline(v = seq(usr[1], usr[2], length = ncol(bigmat) + 1))
abline(h = seq(usr[3], usr[4], length = nrow(bigmat) + 1))

abline(v = seq(usr[1], usr[2], length = ncol(bigmat) + 1)[12:13], col = "purple", lwd = 2)

xx <- seq(usr[1], usr[2], length = ncol(bigmat) + 1)[12:13]
yy <- seq(usr[3], usr[4], length = nrow(bigmat) + 1)[c(1, nrow(bigmat) + 1)]

# LWD IS DEVICE SPECIFIC SO MIGHT HAVE TO ADJUST ON FINAL PDF OUTPUT
YYNUDGE <- 0.002

segments(xx[1], yy[1], xx[2], yy[1], lwd = 2, col = "purple")
segments(xx[1], yy[2]-YYNUDGE, xx[2], yy[2]-YYNUDGE, lwd = 2, col = "purple")

xx <- seq(0, 1, length = ncol(bigmat))
yy <- seq(0, 1, length = nrow(bigmat))

for(cc in 1:length(xx)) {
  for(rr in 1:length(yy)) {
    TeachingDemos::shadowtext(xx[cc], yy[rr], grpmat[rr, cc], col = "grey85", bg = "grey25")
  }
}

addsubplotleg("(a)")



### stag093
n <- which(sapply(info, function(x) x$cur$Tag.ID[1]) == "ZcTag093")

sightmat = info[[n]]$sightmat
focalmat = info[[n]]$focalmat
sightlist = info[[n]]$sightlist
i = info[[n]]$i
xx = info[[n]]$xx
yy = info[[n]]$yy
cur = info[[n]]$cur
rownames_labs = info[[n]]$rownames_labs

rownames_labs <- gsub("Tag", "STag", rownames_labs)

# custom fix for dtag label
rownames_labs[5] <- paste("DTag19_218", rownames_labs[5])

par(mar = c(3.1, 10.1, 3.1, 1.1))
image(t(sightmat + (sightmat & focalmat)), axes = FALSE, col = hcl.colors(10, "blues", rev = TRUE))
axis(2, at = seq(0, 1, length = nrow(sightmat)), lab = rownames_labs, las = 1)

usr <- par('usr')
abline(v = seq(usr[1], usr[2], length = ncol(sightmat) + 1))
abline(h = seq(usr[3], usr[4], length = nrow(sightmat) + 1))

ceecols <- which(colnames(sightmat) %in% sightlist[[2]])
focalrow <- which(rownames(sightmat) == key$Photo.ID[i])

points(xx[ceecols], rep(yy[focalrow], length(ceecols)), pch = 16, cex = 4, col = "grey85")
points(xx[ceecols], rep(yy[focalrow], length(ceecols)), pch = 16, cex = 3.5, col = "purple")

mtext(paste(sub("Zc", "S", cur$Tag.ID), sub("CEE_", "", cur$CEE.Date)), side = 3, line = 2.1)

# calculate +/- days from cee
sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
dates <- as.Date(sight_days, format = "%Y%m%d")
cee_date <- as.Date(as.character(cur$CEE.Date.1), format = "%Y%m%d")
lags <- dates - cee_date
modifiers <- rep("", length(lags))
modifiers[sign(lags) == 1] <- "+"
laglabels <- paste0(modifiers, lags)

axis(1, at = seq(0, 1, length = ncol(sightmat)), lab = laglabels, las = 1)

addsubplotleg("(b)")



### stag097
n <- which(sapply(info, function(x) x$cur$Tag.ID[1]) == "ZcTag097")

sightmat = info[[n]]$sightmat
focalmat = info[[n]]$focalmat
sightlist = info[[n]]$sightlist
i = info[[n]]$i
xx = info[[n]]$xx
yy = info[[n]]$yy
cur = info[[n]]$cur
rownames_labs = info[[n]]$rownames_labs
rownames_labs <- gsub("Tag", "STag", rownames_labs)
rownames_labs <- sub("no5 ", "no. 5 UU", rownames_labs)

# switch 96 to be focal
focalmat[7, ] <- 0
focalmat[3, ] <- 1

par(mar = c(3.1, 10.1, 3.1, 1.1))
image(t(sightmat + (sightmat & focalmat)), axes = FALSE, col = hcl.colors(10, "blues", rev = TRUE))
axis(2, at = seq(0, 1, length = nrow(sightmat)), lab = rownames_labs, las = 1)

usr <- par('usr')
abline(v = seq(usr[1], usr[2], length = ncol(sightmat) + 1))
abline(h = seq(usr[3], usr[4], length = nrow(sightmat) + 1))

ceecols <- which(colnames(sightmat) %in% sightlist[[2]])
focalrow <- which(rownames(sightmat) == "Zca_092_HAT") # hand fix to 96

points(xx[ceecols], rep(yy[focalrow], length(ceecols)), pch = 16, cex = 4, col = "grey85")
points(xx[ceecols], rep(yy[focalrow], length(ceecols)), pch = 16, cex = 3.5, col = "purple")

mtext(paste(sub("Zc", "S", "ZcTag096"), sub("CEE_", "", cur$CEE.Date)), side = 3, line = 2.1)

# calculate +/- days from cee
sight_days <- sapply(strsplit(colnames(sightmat), "_"), '[[', 1)
dates <- as.Date(sight_days, format = "%Y%m%d")
cee_date <- as.Date(as.character(cur$CEE.Date.1), format = "%Y%m%d")
lags <- dates - cee_date
modifiers <- rep("", length(lags))
modifiers[sign(lags) == 1] <- "+"
laglabels <- paste0(modifiers, lags)

axis(1, at = seq(0, 1, length = ncol(sightmat)), lab = laglabels, las = 1)

addsubplotleg("(c)")

dev.off()
