###
# run_gam_plots.r

## some general notes:

# - AIC with a change of at least 2 for selects best model

# - dives are over 33 minutes and over 800 meters and surf(acings) 
#   are less than these numbers (this aligns with our _a prior_ biological
#   definitions of deep diving and bounce diving behavior)


## notes on dive differences in behavior and series:

# - for behavior data *only* if a dive is under 33 minutes and over 800 meters
#   it is not detected as a dive. time-series doesn't have this problem.

# - we can align the time-series so that it also fails to detect these 
#   dives if we want.

# - this is the smallest category of dives so it may be OK to simply 
#   acknowledge the difference.


## filtering gaps:

# - In general, gaps don't matter, but if there isn't 6 hours of continuous
#   data somewhere in the cee_st + 24 hour window, then that cee is excluded.

# - There are a couple of pretty short overall records that are *not* currently
#   filtered out.


## AR:

# - checking is done on the model without exp variables (dive0, surf0)
# - currently gaps are ignored.


## notes on particular tags:

# - A couple of tags, 85 for example, have a bad qq plot, but there isn't
#   much to do about this other than interpret carefully.

# - 60 has a convergence problem.

## logistic notes:

# - sorry the code is a little unnecessarily convoluted at places
#   some of this is regular technical debt, some is because it is designed to
#   deal with arbitrary numbers of covariates for instance.



### helper function for plotting pretty (a), (b), etc. labels
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



### run the models
# AIC criteria
BYTWO <- TRUE

# smooth basis dimension
k <- 3

# model family to fit
model_fam <-  Gamma(link = "log")
# model_fam <- mgcv::tw(link = "log") # an alternative, but fitting isn't good

# load data
dat <- readRDS("../../01_shared_data_products/prep_gam/gam_prepared_dive_dat.rds")
names(dat) <- sapply(dat, function(x) x$param$tag)

# not enough data-- crashes
dat[[24]] <- NULL

# a bunch of clunky lists to hold things
s_mods <- vector(mode = "list", length = length(dat))
bs_mods <- vector(mode = "list", length = length(dat))
d_mods <- vector(mode = "list", length = length(dat))
ds_mods <- vector(mode = "list", length = length(dat))

names(s_mods) <- names(dat)
names(d_mods) <- names(dat)
names(bs_mods) <- names(dat)
names(ds_mods) <- names(dat)

ddexpdur <- vector(mode = "list", length(dat))
dat_processed <- vector(mode = "list", length(dat))

names(ddexpdur) <- names(dat)
names(dat_processed) <- names(dat)

# mega loop
n_tags <- length(dat)
for(i in 1:n_tags) {
  # get the data and smush together the gaps
  input <- dat[[i]]
  divedat <- do.call('rbind', input$dat)
  
  # convert any 0 surface times to 1
  divedat$surface[divedat$surface == 0] <- 1
  
  # get rid of real ship
  divedat <- divedat[divedat$real_ind == 0, ]

  # check to see if we have enough data
  S <- divedat$control_ind + divedat$simulated_ind
  S[S > 0] <- 1
  stretchid <- c(0, cumsum(abs(diff(S)) > 0)) + 1
  stretchid[S == 0] <- 0

  ddexp <- split(divedat, stretchid)[-1]
  ddexpdur[[i]] <- sapply(ddexp, function(x) {
    out <- (max(x$End_Cycle) - min(x$Start)) / 60 / 60
    
    if(nrow(x) > 1) {
      stretchid_gaps <- c(1, cumsum(x$Start[2:nrow(x)] - x$End_Cycle[1:(nrow(x)-1)] > 300) + 1)
      out <- max(sapply(split(x, stretchid_gaps), function(z) {
        (max(z$End_Cycle) - min(z$Start)) / 60 / 60
      }))
    }
 
    out
  })
  
  if(any(ddexpdur[[i]] < 6)) {
    divedat <- divedat[-which(stretchid %in% as.numeric(names(ddexpdur[[i]])[which(ddexpdur[[i]] < 6)])), ]
  }
  
  dat_processed[[i]] <- divedat
  
  # make the exposure variables ordered factors
  sim_ind = factor(divedat$simulated_ind, ordered = TRUE)
  con_ind = factor(divedat$control_ind, ordered = TRUE)

  # set up the variables and labels
  ind <- list(
    SSS = sim_ind,
    SSC = con_ind
  )
  
  ind_var_name <- list(
    SSS = "sim_ind",
    SSC = "con_ind"
  )
  
  ind_var_name_numeric <- list(
    SSS = "simulated_ind",
    SSC = "control_ind"
  )
  
  kmax <- list(
    SSS = min(k, length(unique(divedat$SSS)) - 1),
    SSR = min(k, length(unique(divedat$SSR)) - 1),
    SSC = min(k, length(unique(divedat$SSC)) - 1)
  )
  
  exps <- names(which(sapply(ind, function(x) any(x == 1))))
  
  surf0 <- mgcv::bam(surface ~ 
    s(time, bs = 'ts') +
    s(tsun, bs = 'cc'),
    family = model_fam,
    data = divedat,
    select = TRUE,
    method = "fREML",
    discrete = TRUE
  )
  
  # Account for possible autocorrelation in residuals
  resid_surf0 <- residuals(surf0, type = "response")
  rho_surf0 <- cor(resid_surf0[-length(resid_surf0)], resid_surf0[-1])
  if(rho_surf0 > 0) {
    surf0 <- mgcv::bam(surface ~ 
      s(time, bs = 'ts') +
      s(tsun, bs = 'cc'),
      family = model_fam,
      data = divedat,
      select = TRUE,
      method = "fREML",
      discrete = TRUE,
      rho = rho_surf0
    )
  }
  
  dive0 <- mgcv::bam(dive ~
    s(time, bs = 'ts') +
    s(tsun, bs = 'cc'),
    family = model_fam,
    data = divedat,
    select = TRUE,
    method = "fREML",
    discrete = TRUE
  )
  
  # Account for possible autocorrelation in residuals
  resid_dive0 <- residuals(dive0, type = "response")
  rho_dive0 <- cor(resid_dive0[-length(resid_dive0)], resid_dive0[-1])
  if(rho_dive0 > 0) {
    dive0 <- mgcv::bam(dive ~ 
      s(time, bs = 'ts') +
      s(tsun, bs = 'cc'),
      family = model_fam,
      data = divedat,
      select = TRUE,
      method = "fREML",
      discrete = TRUE,
      rho = rho_dive0
    )
  }
  
  # add terms in for exposure
  s_mods_tmp <- list()
  d_mods_tmp <- list()
  
  bs_mods[[i]] <- surf0
  ds_mods[[i]] <- dive0
  
  s_mods_tmp[[1]] <- surf0
  d_mods_tmp[[1]] <- dive0
  
  terms <- vector(mode = "list", length = length(exps))
  count <- 1
  
  if(length(exps) > 0) {
    for(n in 1:length(exps)) {
      if(kmax[[exps[n]]] < 3) {
        terms[[count]] <- ind_var_name[[exps[n]]]
      } else {
        terms[[count]] <- paste0(
          "s(", exps[n], ", ",
          "by = ", ind_var_name_numeric[[exps[n]]],
          ", bs = 'ts', k = ",
          kmax[[exps[n]]],
          ", m = 1, pc = 1440)"
        )
      }
      
      count <- count + 1
    }
  
    # possible combinations
    combos <- list()
    
    for(n in 1:length(terms)) {
      combos[[n]] <- combn(length(terms), n)
    }
    
    modupdate <- list()
    count <- 1
    for(n in 1:length(combos)) {
      curcombo <- combos[[n]]
      
      for(b in 1:ncol(curcombo)) {
        modupdate[[count]] <- paste(terms[curcombo[, b]], collapse = " + ")
        count <- count + 1
      }
    }
    
    # do it
    s_count <- length(s_mods_tmp) + 1
    d_count <- length(d_mods_tmp) + 1
    for(n in 1:length(modupdate)) {
      s_try <- tryCatch(update(s_mods_tmp[[1]], paste("~. +", modupdate[[n]])), error = function(e) warning(e))
      d_try <- tryCatch(update(d_mods_tmp[[1]], paste("~. +", modupdate[[n]])), error = function(e) warning(e))
      
      if("gam" %in% class(s_try)) {
        s_mods_tmp[[s_count]] <- s_try
        s_count <- s_count + 1
      }
      
      if("gam" %in% class(d_try)) {
        d_mods_tmp[[d_count]] <- d_try
        d_count <- d_count + 1
      }
    }
  }  
  # apply AIC rule

  # surf
  aics <- sapply(s_mods_tmp, AIC)
  s_dis <- 1
    
  if(BYTWO) {
    dese <- aics[1] - aics >= 2
    if(any(dese)) {
      minaic <- min(aics[dese])
      s_dis <- which(aics == minaic)
    }
  } else {
    # subtract a tiny amount from AIC for simplest model to
    # ensure it is selected if AICs are equal
    aics[1] <- aics[1] - 0.01
    s_dis <- which(aics == min(aics))
  }
  
  # dive 
  aics <- sapply(d_mods_tmp, AIC)
  d_dis <- 1
    
  if(BYTWO) {
    dese <- aics[1] - aics >= 2
    if(any(dese)) {
      minaic <- min(aics[dese])
      d_dis <- which(aics == minaic)
    }
  } else {
    # subtract a tiny amount from AIC for simplest model to
    # ensure it is selected if AICs are equal
    aics[1] <- aics[1] - 0.01
    d_dis <- which(aics == min(aics))
  }
  
  s_mods[[i]] <- s_mods_tmp[[s_dis]]
  d_mods[[i]] <- d_mods_tmp[[d_dis]]
  
}

# make significance table
sigdat <- data.frame(
  id = names(dat), 
  surf_sim = NA, 
  dive_sim = NA, 
  surf_con = NA, 
  dive_con = NA
)

# extract sig function
extract_sig <- function(x, cee) {
  out <- NA
  
  if("gam" %in% class(x)) {
    ff <- x$formula
    out <- length(grep(cee, ff)) > 0
  }
  
  if(length(out) == 0) out <- NA
  
  out
}

# extract p values
sigdat$surf_sim <- sapply(s_mods, extract_sig, "sim")
sigdat$dive_sim <- sapply(d_mods, extract_sig, "sim")
sigdat$surf_con <- sapply(s_mods, extract_sig, "con")
sigdat$dive_con <- sapply(d_mods, extract_sig, "con")

# figuring out again which cees are actually relevant
# load additional metadata
meta <- read.table("../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
ceemeta <- read.table("../../00_data_input/cee/cee_metadata_flat.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# clean up dates
ceemeta$start_datenum <- as.numeric(as.POSIXct(paste(ceemeta$date_YYYYMMDD, ceemeta$start_HHMMSS), format = "%Y%m%d %H%M%S", tz = 'utc'))
st <- ceemeta$start_datenum
en <- st + 1* 24 * 60 * 60

# find the considered cees
considered_cees <- vector()
considered_dids <- vector()

for(i in 1:length(dat_processed)) {
  xx <- dat_processed[[i]]
  res <- vector(length = length(st))
  for(n in 1:length(st)) {
    res[n] <- any(findInterval(xx$Start, c(st[n], en[n])) == 1) 
  }
  
  considered_cees <- c(considered_cees, ceemeta$cee_id[res])
  considered_dids <- c(considered_dids, rep(names(dat_processed)[i], length(which(res))))
}

considered <- data.frame(deployid = considered_dids, cee_id = considered_cees)
considered$type <- ceemeta$cee_type[match(considered_cees, ceemeta$cee_id)]

key <- paste(considered$deployid, considered$cee_id)
tab <- paste(meta$deployid, meta$cee_id)

considered$spl <- meta$median_spl[match(key, tab)]
considered$rng <- meta$median_range_km[match(key, tab)]
considered$spl_bin <- meta$rl_bin_spl[match(key, tab)]



### merge considered with the results and output a table

# set a seed for the jitter
set.seed(0000000311828578)

sigdat_mer <- merge(sigdat, considered, by.x = "id", by.y = "deployid")

# clean up
sigdat_mer$surf_con[sigdat_mer$type != "control"] <- NA
sigdat_mer$dive_con[sigdat_mer$type != "control"] <- NA
sigdat_mer$surf_sim[sigdat_mer$type != "simulated mfas"] <- NA
sigdat_mer$dive_sim[sigdat_mer$type != "simulated mfas"] <- NA

# output data
write.table(sigdat_mer, "gam_results_aic.csv", sep = ',', row.names = FALSE)



### make the plots
pdf("gam_summary.pdf", width = 8, height = 6)
par(mfrow = c(2, 3), family = "serif")
s93 <- s_mods[['ZcTag093_DUML']]
d93 <- d_mods[['ZcTag093_DUML']]
e93 <- dat_processed[['ZcTag093_DUML']]
bs93 <- bs_mods[['ZcTag093_DUML']]
bd93 <- ds_mods[['ZcTag093_DUML']]

nsim <- 10000



### subplot a
newt <- seq(min(e93$time), max(e93$time), length = 1000)

newd <- data.frame(
  time = newt,
  tsun = approx(e93$time, e93$tsun, newt)$y,
  SSS = approx(e93$time, e93$SSS, newt)$y,
  simulated_ind = approx(e93$time, e93$simulated_ind, newt)$y
)

br <- MASS::mvrnorm(nsim, coef(s93), s93$Vp)
Xp <- predict(s93, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

ylims <- c(min(min(cl2), 0), max(max(cl2), max(e93$surface)))

plot(e93$time[e93$simulated_ind != 1]/60, e93$surface[e93$simulated_ind != 1], 
  col = "grey", pch = 16, 
  ylim = ylims, 
  xlab = 'running time (hours)', ylab = "IDDI (min)", 
  las = 1
)

points(e93$time[e93$simulated_ind == 1]/60, e93$surface[e93$simulated_ind == 1],
  col = 2, pch = 16
)

lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 2, lwd = 2)

# add to plot
lines(newd$time/60, cl2[1, ], col = 2, lty = 2)
lines(newd$time/60, cl2[2, ], col = 2, lty = 2)

# do it again with no exposure term
br <- MASS::mvrnorm(nsim, coef(bs93), bs93$Vp)
Xp <- predict(bs93, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 1, lwd = 1)

# add to plot
lines(newd$time/60, cl2[1, ], col = 1, lty = 3)
lines(newd$time/60, cl2[2, ], col = 1, lty = 3)

title("STag093 19_03")
legend("topleft", legend = "change detected", bty = 'n', text.col = "purple")

addsubplotleg("(a)")



### subplot b
newt <- seq(min(e93$time), max(e93$time), length = 1000)

newd <- data.frame(
  time = newt,
  tsun = approx(e93$time, e93$tsun, newt)$y,
  SSS = approx(e93$time, e93$SSS, newt)$y,
  simulated_ind = approx(e93$time, e93$simulated_ind, newt)$y
)

newd$simulated_ind[newd$simulated_ind > 0] <- 1

br <- MASS::mvrnorm(nsim, coef(d93), d93$Vp)
Xp <- predict(d93, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

ylims <- c(min(min(cl2), 33), max(max(cl2), max(e93$dive)))

plot(e93$time[e93$simulated_ind != 1]/60, e93$dive[e93$simulated_ind != 1], col = "grey", pch = 16, ylim = ylims, ylab = "deep-dive duration (min)", xlab = "running time (hours)", las = 1)
points(e93$time[e93$simulated_ind == 1]/60, e93$dive[e93$simulated_ind == 1], col = 2, pch = 16)

# do it again with no exposure term
br <- MASS::mvrnorm(nsim, coef(bd93), bd93$Vp)
Xp <- predict(bd93, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 1, lwd = 1)

# add to plot
lines(newd$time/60, cl2[1, ], col = 1, lty = 3)
lines(newd$time/60, cl2[2, ], col = 1, lty = 3)

title("STag093 19_03")
addsubplotleg("(b)")

# add a legend for the model plots
plot(1, 1, type = 'n', xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(
  "obs.", 
  "obs. during exposure", 
  "pred baseline model", 
  "pred exposure model", 
  "95% CI baseline model", 
  "95% CI exposure model"),
  pch = c(16, 16, NA, NA, NA, NA),
  col = c("grey", 2, 1, 2, 1, 2),
  lty = c(NA, NA, 1, 1, 3, 2)
)



### subplot c
s96 <- s_mods[['ZcTag096_DUML']]
d96 <- d_mods[['ZcTag096_DUML']]
e96 <- dat_processed[['ZcTag096_DUML']]
bs96 <- bs_mods[['ZcTag096_DUML']]
bd96 <- ds_mods[['ZcTag096_DUML']]

nsim <- 10000

newt <- seq(min(e96$time), max(e96$time), length = 1000)

newd <- data.frame(
  time = newt,
  tsun = approx(e96$time, e96$tsun, newt)$y,
  SSS = approx(e96$time, e96$SSS, newt)$y,
  simulated_ind = approx(e96$time, e96$simulated_ind, newt)$y
)

br <- MASS::mvrnorm(nsim, coef(s96), s96$Vp)
Xp <- predict(s96, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

ylims <- c(min(min(cl2), 0), max(max(cl2), max(e96$surface)))

plot(e96$time[e96$simulated_ind != 1]/60, e96$surface[e96$simulated_ind != 1], col = "grey", pch = 16, ylim = ylims, ylab = "IDDI (min)", xlab = "running time (hours)", las = 1)
points(e96$time[e96$simulated_ind == 1]/60, e96$surface[e96$simulated_ind == 1], col = 2, pch = 16)

# do it again with no exposure term
br <- MASS::mvrnorm(nsim, coef(bs96), bs96$Vp)
Xp <- predict(bs96, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 1, lwd = 1)

# add to plot
lines(newd$time/60, cl2[1, ], col = 1, lty = 3)
lines(newd$time/60, cl2[2, ], col = 1, lty = 3)

title("STag096 19_04")
addsubplotleg("(c)")



### subplot (d)

newt <- seq(min(e96$time), max(e96$time), length = 1000)

newd <- data.frame(
  time = newt,
  tsun = approx(e96$time, e96$tsun, newt)$y,
  SSS = approx(e96$time, e96$SSS, newt)$y,
  simulated_ind = approx(e96$time, e96$simulated_ind, newt)$y
)

newd$simulated_ind[newd$simulated_ind != 1] <- 0

br <- MASS::mvrnorm(nsim, coef(d96), d96$Vp)
Xp <- predict(d96, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

ylims <- c(min(min(cl2), 33), max(max(cl2), max(e96$dive)))

plot(e96$time[e96$simulated_ind != 1]/60, e96$dive[e96$simulated_ind != 1], col = "grey", pch = 16, ylim = ylims, ylab = "deep-dive duration (min)", xlab = "running time (hours)", las = 1)
points(e96$time[e96$simulated_ind == 1]/60, e96$dive[e96$simulated_ind == 1], col = 2, pch = 16)
lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 2, lwd = 2)

# add to plot
lines(newd$time/60, cl2[1, ], col = 2, lty = 2)
lines(newd$time/60, cl2[2, ], col = 2, lty = 2)

# do it again with no exposure term
br <- MASS::mvrnorm(nsim, coef(bd96), bd96$Vp)
Xp <- predict(bd96, newd, type = 'lpmatrix')
res <- exp(Xp %*% t(br))
cl2 <- apply(res, 1, quantile, c(0.025, 0.975))

lines(newd$time/60, apply(res, 1, quantile, 0.5), col = 1, lwd = 1)

# add to plot
lines(newd$time/60, cl2[1, ], col = 1, lty = 3)
lines(newd$time/60, cl2[2, ], col = 1, lty = 3)

title("STag096 19_04")
legend("topleft", legend = "change detected", bty = 'n', text.col = "purple")
addsubplotleg("(d)")



### subplot e
# binned summary
sigdat_mer <- sigdat_mer[!is.na(sigdat_mer$spl_bin), ]

# melt
sigdat_mer_tall <- reshape2::melt(sigdat_mer,
  id.vars = c('id', 'cee_id', 'spl', 'rng', 'spl_bin', 'type'),
  variable.name = 'test',
  value.name = 'sig'
)

# remove unwanted combinations
keep1 <- grepl("con", sigdat_mer_tall$test) & sigdat_mer_tall$type == "control"
keep2 <- grepl("sim", sigdat_mer_tall$test) & sigdat_mer_tall$type == "simulated mfas"

sigdat_mer_tall <- sigdat_mer_tall[keep1 | keep2, ]

# this taken care of by type
sigdat_mer_tall$test <- sapply(strsplit(as.character(sigdat_mer_tall$test), "_"), '[[', 1)

# order the levels for nicer plot
sigdat_mer_tall$spl_bin <- factor(sigdat_mer_tall$spl_bin, levels = c('control', '<100', '100-120', '>120'))

ddtot <- tapply(
  sigdat_mer_tall$sig[sigdat_mer_tall$test == "dive"],
  sigdat_mer_tall$spl_bin[sigdat_mer_tall$test == "dive"],
  function(x) length(which(x))
)

sstot <- tapply(
  sigdat_mer_tall$sig[sigdat_mer_tall$test == "surf"],
  sigdat_mer_tall$spl_bin[sigdat_mer_tall$test == "surf"],
  function(x) length(which(x))
)

totot <- tapply(
  sigdat_mer_tall$sig,
  sigdat_mer_tall$spl_bin,
  function(x) length(which(!is.na(x)))
)

yyy <- rbind(ddtot/totot, sstot/totot)

yydive <- yyy[1, ]
yysurf <- yyy[2, ]
xx1 <- 1:4
xx2 <- seq(1.25, length = 4, by = 1)

plot(0, 0, axes = FALSE, type = 'n', xlab = "RL bin (dB SPL)", ylab = "percentage response identified", xlim = c(0.75, 4.75), ylim = c(0, 1), main = "summary of all tags")
box()
points(xx1, yydive, pch = 16, col = 3)
points(xx2, yysurf, pch = 16, col = 1)
axis(1, at = apply(rbind(xx1, xx2), 2, mean), labels = names(yydive))
axis(2, las = 1)
abline(h = pretty(c(0, 1)), lty = 3, col = "grey75")
legend("topleft", legend = c("deep-dive duration", "IDDI"), pch = 16, col = c(3, 1), bty = 'n')
text(apply(rbind(xx1, xx2), 2, mean), rep(.45, 5), paste("n =", totot/2))
# axis(1, at = apply(rbind(xx1, xx2), 2, mean), line = 2.1, labels = paste("n =", totot/2), lty = 0, )
addsubplotleg("(e)")

dev.off()



### results table
sigdat_mer

# double check
length(which(!is.na(sigdat_mer$surf_con)))
length(which(!is.na(sigdat_mer$dive_con)))

length(which(!is.na(sigdat_mer$dive_sim[sigdat_mer$spl_bin == ">120"])))
length(which(!is.na(sigdat_mer$surf_sim[sigdat_mer$spl_bin == ">120"])))

length(which(!is.na(sigdat_mer$dive_sim[sigdat_mer$spl_bin == "100-120"])))
length(which(!is.na(sigdat_mer$surf_sim[sigdat_mer$spl_bin == "100-120"])))

length(which(!is.na(sigdat_mer$dive_sim[sigdat_mer$spl_bin == "<100"])))
length(which(!is.na(sigdat_mer$surf_sim[sigdat_mer$spl_bin == "<100"])))
