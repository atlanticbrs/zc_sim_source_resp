###
# prep prh

# libraries
library(R.matlab)
library(signal)
source("find_dives.r")

# metadata
dtagontime <- as.numeric(as.POSIXct("2019-08-06 15:38:20 UTC", tz = "UTC"))
tagoff_cutoff <- as.numeric(as.POSIXct("2019-08-06 20:16:40 UTC", tz = "UTC"))

# load prh
prh <- R.matlab::readMat("zc19_218aprh.mat")

# calculate msa
prh$msa <- 9.81 * matrix(abs(sqrt(rowSums(prh$A^2)) - 1),
  nrow(prh$A), 1
)
# prh$msa <- abs(sqrt(prh$A^2 %*% cbind(rep(1, 3))) - 1) # faster than rowSUMS???
# yes rowSums is slower, it does a bunch of checking and has more features

# calculate a time vector
# prh$fs is actually a 1x1 array and R has depreciated recyucling it (incr/prh$fs)
# as.vector fixes the warning
prh$t_cues <- matrix(1:nrow(prh$p) / as.vector(prh$fs),
  nrow(prh$p), 1
)

prh$t <- matrix(dtagontime + prh$t_cues,
  nrow(prh$t_cues), 1
)

# filter to tagoff time
dese <- prh$t <= tagoff_cutoff

for(i in 1:length(prh)) {
  if(length(prh[[i]]) != 1) {
    prh[[i]] <- prh[[i]][dese, , drop = FALSE]
  }
}

# 
# windows()
# par(mfrow = c(2, 1))
# plot(1:nrow(prh$p), -prh$p, type = 'l')
# plot(1:nrow(prh$p), msa*9.81, type = 'l', ylim = c(0, 10), lwd = 3)

dive_cues <- find_dives(prh$p, 50, as.vector(prh$fs), surface = 1, findall = FALSE)

# translate cues to indices
dive_cues$st <- match(dive_cues$start, prh$t_cues)
dive_cues$en <- match(dive_cues$end, prh$t_cues)

prh$dive_cues <- dive_cues

saveRDS(prh, file = "Zc19_218a_prh.rds")
