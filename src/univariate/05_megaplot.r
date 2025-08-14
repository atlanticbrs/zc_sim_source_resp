###
# megaplot
# make a multipanel plot summarizing univariate data



### constants
NM_PER_KM <- 1.852001
MAX_MODELLING_DISTANCE_RL_NM <- 70
MAX_MODELLING_DISTANCE_RL_KM <- MAX_MODELLING_DISTANCE_RL_NM * NM_PER_KM



### required libraries and helper functions
library(corrplot)

source("01_helperfunctions/findgaps.r")
source("01_helperfunctions/plot_series.r")

# helper to make (a), (b), ..., etc. legends on par multipanel plots
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

# histo plot function
make_metric_histo <- function(cur, ...) {
  obs <- cur$obs
  bas <- cur$bas
  
  # make xlim to include the exp dive
  # xlims <- range(c(obs, bas), na.rm = TRUE) * c(0.8, 1.2) # 20% slop
  
  hh <- hist(bas,
             # xlim = xlims,
             main = "",
             xlab = "IDDI (minutes)",
             ylab = "frequency",
             yaxt = 'n',
             ...
  )
  
  axis(2, las = 1)
  
  # make obs arrow
  arrowheight <- max(hh$counts) * 0.40
  arrows(obs, arrowheight, obs, 0, length = 0.15, col = "purple", lwd = 2)
  
  # make labels
  bas_nona <- bas[!is.na(bas)]
  perc = round(length(which(obs > bas_nona)) / length(bas_nona), 3)
  # legend("topright", c(paste("CEE", cur$info$ceeid), cur$info$deployid, cur$info$label_selected, cur$info$settings, paste("perc =", perc)))
  
  box()
}



### make plot
pdf("uniplot.pdf", width = 10, height = 8)
# png("uniplot.png", width = 10*300, height = 8*300, res = 300)

lay <- matrix(c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7), 4, 3, byrow = FALSE)
layout(lay)
par(mar = c(4.1, 5.1, 2.1, 1.1), family = "serif", xpd = FALSE)



### dive plots
plot_inputs <- readRDS("02_intermediate_outputs/plot_inputs.rds")
plot_type <- sapply(plot_inputs, '[[', 'what')
deployid <- sapply(plot_inputs, function(x) x$info$deployid)

curinput <- plot_inputs[[which(plot_type == "IDDI" & deployid == "ZcTag093_DUML")]]
cee93 <- curinput$info

curinput <- plot_inputs[[which(plot_type == "IDDI" & deployid == "ZcTag096_DUML")]]
cee96 <- curinput$info

# 93 and 89 dive profiles
# load data
s93 <- read.table("../../01_shared_data_products/filter_sattag/series/ZcTag093_DUML_180763_L2_gonio-pressuresensor-Series.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
s96 <- read.table("../../01_shared_data_products/filter_sattag/series/ZcTag096_DUML_180774_L2_gonio-pressuresensor-Series.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# 93
stx <- cee93$ceest - 12 * 60 * 60
enx <- stx + 24 * 60 * 60

plot_series(s93, las = 1, xaxt = 'n', yaxt = 'n', xlim = c(stx, enx), xlab = expression(t - t[exposure_init]~(hours)), show_gaps = FALSE, ylab = "depth (meters)")

xx <- seq(stx, enx, by = 60*60)
xxlab <- (xx - cee93$ceest) / 60 / 60
axis(1, at = xx, lab = xxlab)
axis(2, at = c(0, -500, -1000, -1500), lab = c(0, 500, 1000, 1500), las = 1)

rect(cee93$ceest, -10000, cee93$ceeen, 10000, col = rgb(0.25, 0, 1, .25), border = NA)

title("STag093 19_03")
addsubplotleg("(a)")


# 96
stx <- cee96$ceest - 12 * 60 * 60
enx <- stx + 24 * 60 * 60

plot_series(s96, las = 1, xaxt = 'n', yaxt = 'n', xlim = c(stx, enx), xlab = expression(t - t[exposure_init]~(hours)), show_gaps = FALSE, ylab = "depth (meters)")

xx <- seq(stx, enx, by = 60*60)
xxlab <- (xx - cee96$ceest) / 60 / 60
axis(1, at = xx, lab = xxlab)

axis(2, at = c(0, -500, -1000, -1500), lab = c(0, 500, 1000, 1500), las = 1)

rect(cee96$ceest, -10000, cee96$ceeen, 10000, col = rgb(0.25, 0, 1, .25), border = NA)

title("STag096 19_04")
addsubplotleg("(b)")



### histograms
plot_type <- sapply(plot_inputs, '[[', 'what')
deployid  <- sapply(plot_inputs, function(x) x$info$deployid)
ceeid <- sapply(plot_inputs, function(x) x$info$ceeid)
ceeyear <- as.numeric(sapply(strsplit(ceeid, '_'), '[[', 1)) + 2000
ceetype <- sapply(plot_inputs, function(x) x$info$ceetype)

# plot IDDI for two tags
curinput <- plot_inputs[[which(plot_type == "IDDI" & deployid == "ZcTag093_DUML")]]
cee93 <- curinput$info
make_metric_histo(curinput, nclass= 15, xlim = c(0, 425))
title("STag093 19_03")
text(curinput$obs, 10, "q = 1.00", col = "purple")
addsubplotleg("(c)")

curinput <- plot_inputs[[which(plot_type == "IDDI" & deployid == "ZcTag096_DUML")]]
cee96 <- curinput$info
make_metric_histo(curinput, nclass = 20, xlim = c(0, 425))
title("STag096 19_04")
text(curinput$obs, 10, "q = 0.99", col = "purple")
addsubplotleg("(d)")



### boxplot
# load rls and make a key
rls <- read.table("../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# filter and make key
rls <- rls[!is.na(rls$rl_bin_spl), ]
rlkey_rls <- paste(rls$deployid, rls$cee_id)

exp_type <- sapply(plot_inputs, function(x) x$info$ceetype)
dese <- which(exp_type %in% c("simulated mfas", "control"))

dat <- plot_inputs[dese]
datperc <- sapply(dat, function(x) {
  bas <- x$bas
  bas_nona <- bas[!is.na(bas)]
  length(which(x$obs > bas_nona)) / length(bas_nona)
})

sig <- datperc > 0.975 | datperc < 0.025
byanimalbycee <- sapply(dat, function(x) paste(x$info$deployid, x$info$ceeid))

sig_split <- split(sig, byanimalbycee)
percsigmetric <- sapply(sig_split, function(x) length(which(x)) / length(x))
nmetrics <- sapply(sig_split, function(x) length(x))

rlkey_dat <- names(percsigmetric)

rlbins <- rls$rl_bin_spl[match(rlkey_dat, rlkey_rls)]
rlbins <- factor(rlbins, levels = c("control", "<100", "100-120", ">120"))

# calculate medians
meds <- as.data.frame(t(rbind(
  med = tapply(percsigmetric, rlbins, median),
  low = tapply(percsigmetric, rlbins, quantile, 0.25),
  hih = tapply(percsigmetric, rlbins, quantile, 0.75)
)))

# calculate x
xx <- as.numeric(rlbins)

plot(jitter(xx), percsigmetric, axes = FALSE, ylab = "percentage of metrics identified", xlab = "RL bins (dB SPL)", col = rgb(0, 0, 1, .25), pch = 16, cex = log(nmetrics))
axis(1, at = 1:nrow(meds), lab = rownames(meds))
axis(2, las = 1)
box()

points(1:nrow(meds), meds$med, pch = 16)
segments(1:nrow(meds), meds$low, 1:nrow(meds), meds$hih)

addsubplotleg("(e)")



### bar plot
deployid <- vector()
cee <- vector()
cee_type <- vector()
perc <- vector()
sig <- vector()
what <- vector()
settings <- vector()

for(i in 1:length(plot_inputs)) {
  cur <- plot_inputs[[i]]
  
  # calculate perc and sig
  obs <- cur$obs
  bas <- cur$bas
  
  bas_nona <- bas[!is.na(bas)]
  perc[i] = length(which(obs > bas_nona)) / length(bas_nona)
  
  sig[i] <- perc[i] > 0.975 | perc[i] < 0.025
  
  # add metadata
  deployid[i] <- cur$info$deployid
  cee[i] <- cur$info$ceeid
  cee_type[i] <- cur$info$ceetype
  what[i] <- cur$what
  settings[i] <- cur$info$settings
}

# make data frame
dat <- data.frame(deployid, cee, cee_type, rl_cat = rep(NA, length(deployid)), perc, sig, what, settings)

# merge in rlcat
dat$rl_cat <- rls$rl_bin_spl[match(
  paste(dat$deployid, dat$cee),
  paste(rls$deployid, rls$cee_id)
)]

# check for NAs
dat$deployid[is.na(dat$rl_cat)]

# delete rows with NAs
dat <- dat[!is.na(dat$rl_cat), ]

# calculate yes no in each category
dat2 <- dat[, c('what', 'sig', 'rl_cat')]
ag <- aggregate(sig ~ ., data = dat2, function(x) length(which(x)))
sig_yes <- xtabs(sig ~ ., data = ag)

ag <- aggregate(sig ~ ., data = dat2, function(x) length(which(!x)))
sig_no <- xtabs(sig ~ ., data = ag)

# row order
oo <- order(rowSums(sig_yes))
sig_yes <- sig_yes[oo, ]

# labels
labs <- rownames(sig_yes)
labs <- sub("_", " ", sub("_", " ", labs))

# change bounce to shallower dive
labs <- sub("bounce", "shallower dive", labs)

sig_tmp <- reshape2::melt(sig_yes)
sig_tmp$sig <- "yes"
sig_tmp2 <- reshape2::melt(sig_no)
sig_tmp2$sig <- "no"

sig <- rbind(sig_tmp, sig_tmp2)

# col order
oo <- c(4, 1, 3, 2)

bp <- barplot(sig_yes[, oo], xlab = "RL bins (dB SPL)", ylab = "number of metrics identified", las = 1, col = c(palette.colors(), "white"), ylim = c(0, 30))
box()
legend("topleft", legend = c(gsub("bounce", "shallower dive", gsub("maxdepth", "max. depth", gsub("_", " ", rownames(sig_yes))))), pt.bg = c(palette.colors(), "white"), pch = 22, bty = 'n')

yes <- colSums(sig_yes)[oo]
no <- colSums(sig_no)[oo]

yesnolabs <- paste0("n = ", yes, "/", yes+no)
text(bp, yes + 1, yesnolabs)

addsubplotleg("(f)")



### grid
inputkey <- sapply(plot_inputs, function(x) paste(x$info$deployid, x$info$ceeid))
rlkey <- paste(rls$deployid, rls$cee_id)

# all sim + controls
curinput <- plot_inputs[inputkey %in% rlkey]

# make columns and rows
coltypes <- sort(unique(sapply(curinput, '[[', 'what')))
rowtypes <- sort(unique(sapply(curinput, function(x) paste(x$info$deployid, x$info$ceeid))))

m <- matrix(NA, nrow = length(rowtypes), ncol = length(coltypes))
dimnames(m) <- list(rowtypes, coltypes)

for(i in 1:length(curinput)) {
  cur <- curinput[[i]]
  colid <- cur$what
  rowid <- paste(cur$info$deployid, cur$info$ceeid)
  
  # calc perc
  obs <- cur$obs
  bas <- cur$bas
  bas_nona <- bas[!is.na(bas)]
  perc = round(length(which(obs > bas_nona)) / length(bas_nona), 3)
  
  # assign perc
  m[rowid, colid] <- perc
}

# match to rl bin
bins <- rls$rl_bin_spl[match(rownames(m), rlkey)]
bins <- factor(bins, levels = c("control", "<100", "100-120", ">120"))

# clean up the row labels
rowtypes <- sub("Zc", "S", rowtypes)
rowtypes <- sub("_DUML", "", rowtypes)
rowtypes <- paste0(rowtypes, " (", bins, ")")
coltypes <- gsub("_", " ", coltypes)
coltypes <- gsub("maxdepth", "max. depth", coltypes)
coltypes <- gsub("bounce", "shallower dive", coltypes)

dimnames(m) <- list(rowtypes, coltypes)

oo <- order(rowSums(m > 0.975 | m < 0.025, na.rm = TRUE))
m <- m[oo, ]
bins <- bins[oo]

oo <- order(bins)
m <- m[rev(oo), ]

oo <- order(colSums(m > 0.975 | m < 0.025, na.rm = TRUE))
m <- m[, rev(oo)]

rownames(m)[grep(">120", rownames(m))[1]] <- paste(   ">120 dB SPL     ", rownames(m)[grep(">120", rownames(m))[1]])
rownames(m)[grep("100-120", rownames(m))[1]] <- paste("100-120 dB SPL  ", rownames(m)[grep("100-120", rownames(m))[1]])
rownames(m)[grep("<100", rownames(m))[1]] <- paste(   "<100 dB SPL     ", rownames(m)[grep("<100", rownames(m))[1]])
rownames(m)[grep("control", rownames(m))[1]] <- paste("control         ", rownames(m)[grep("control", rownames(m))[1]])

rownames(m) <- gsub(".\\(100-120\\)", "", rownames(m))
rownames(m) <- gsub(".\\(>120\\)", "", rownames(m))
rownames(m) <- gsub(".\\(<100\\)", "", rownames(m))
rownames(m) <- gsub(".\\(control\\)", "", rownames(m))

corrplot::corrplot(m, method = "color", is.corr = FALSE, col.lim = c(0, 1), col = corrplot::COL2('RdBu', 200), tl.col = "black", tl.srt = 60, na.label = "/", tl.cex = 0.75, cl.pos = 'r', cl.cex = 0.75)

xseq <- seq(0.5, 0.5 + ncol(m), by = 1)
yseq <- seq(0.5, 0.5 + nrow(m), by = 1)

segments(xseq, yseq[1], xseq, yseq[length(yseq)], col = "black", lwd = 2)
segments(xseq[1], yseq, xseq[length(xseq)], yseq, col = "black", lwd = 2)

sigmat <- m > 0.975 | m < 0.025
sigmat_symbol <- sigmat
sigmat_symbol[sigmat]  <- "*"
sigmat_symbol[!sigmat] <- ""

xmat <- matrix(1:ncol(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
ymat <- matrix(nrow(m):1, nrow = nrow(m), ncol = ncol(m), byrow = FALSE)

text(xmat, ymat, sigmat_symbol, col = "grey85", cex = 1.5)
text(xmat[!sigmat], ymat[!sigmat], round(m[!sigmat], 2), col = "black", cex = 0.75)

addsubplotleg("(g)")

dev.off()