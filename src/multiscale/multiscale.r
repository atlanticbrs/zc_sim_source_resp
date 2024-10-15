###
# multiscale.r
# combine datastreams from sattag and dtag



### libraries and helper functions
source("01_helperfunctions/distance.r")

# find the nearest time
matchtimes <- function(t1, t2) {
  findInterval(t1, c(-Inf, head(t2, -1)) + c(0, diff(t2)/2))
}

# pretty subplot (a), (b), ..., etc. legends for par multipanel plots
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



### constants
dtagontime <- as.numeric(as.POSIXct("2019-08-06 15:38:20 UTC", tz = "UTC"))
tagoff_cutoff <- 1565122600



### load data
stag <- read.table("00_data/ZcTag093_DUML_series.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
dtag <- readRDS("00_data/prh_shearer/Zc19_218a_prh.rds")
meta <- read.table("../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
load("00_data/tracks.rdata")
track_datenum <- as.numeric(track_datetime)



### make plot
# filter the track to the dive data
xx <- xx[, track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]
yy <- yy[, track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]
xx_m <- xx_m[track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]
yy_m <- yy_m[track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]
track_dese <- track_dese[track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]
track_datenum <- track_datenum[track_datenum >= min(stag$Date) & track_datenum <= max(stag$Date)]


# sort stag
stag <- stag[order(stag$Date), ]

# make a window around CEE
ceest <- meta$cee_st[meta$cee_id == "19_03"][2]
winbacklag_hr <- 12
winfrontlag_hr <- 12

win <- c(ceest - winbacklag_hr * 60 * 60, ceest + winfrontlag_hr * 60 * 60)

# put in gaps in the stag
stag_period <- diff(stag$Date[1:2])
stag$stretchid <- c(0, cumsum(diff(stag$Date) > stag_period)) + 1

stag_filt <- stag[stag$Date >= win[1] & stag$Date <= win[2], ]

lay <- matrix(c(
  1, 1, 1, 1,
  2, 2, 2, 2,
  3, 3, 3, 3,
  3, 3, 3, 3,
  4, 4, 4, 4,
  4, 4, 4, 4,
  4, 4, 4, 4,
  5, 5, 5, 5,
  5, 5, 5, 5,
  6, 6, 6, 6,
  6, 6, 6, 6,
  6, 6, 6, 6,
  7, 7, 7, 7,
  7, 7, 7, 7
), 14, 4, byrow = TRUE)

# png("multiscale.png", width = 10 * 320, height = 10 * 320, res = 320)
pdf("multiscale.pdf", width = 10, height = 10)

layout(lay)
par(mar = c(0.1, 11.1, 0.1, 0.1), family = "serif")

# track filt
track_filt <- track_datenum >= win[1] & track_datenum <= win[2]

xlim_prime <- min(ceest - min(stag$Date), max(stag$Date) - ceest)
xlims <- c(ceest - xlim_prime, ceest + xlim_prime)

ylims <- range(unlist(yy))

# lat
plot(track_datenum, yy_m, type = 'n', axes = FALSE, xlab = '', ylab = '', ylim = ylims, xlim = xlims)
abline(h = c(35.40, 35.5, 35.6), col = "grey85", lty = 2)

for(i in 1:100) {
  lines(track_datenum, yy[i, ], col = rgb(0, 0, 0, .02))
}

points(track_datenum[!track_filt & track_dese], yy_m[!track_filt & track_dese], pch = 16, cex = .5)
# axis(2, las = 1)
axis(2, las = 1, at = c(35.40, 35.50, 35.6), labels = c("35\u00B0 24\'", "35\u00B0 30\'", "35\u00B0 36\'"))
# legend("topleft", legend = "latitude", bty = 'n')
box(lwd = 1.5)
# rect(win[1], -10000, win[2], 10000, col = rgb(0, 0, 1, .25), border = NA)
abline(v = ceest, col = "#EF476F", lty = 2, lwd = 2)

lines(track_datenum[track_filt & track_dese], yy_m[track_filt & track_dese], col = "#6186DC")
# points(track_datenum[track_filt & track_dese], yy_m[track_filt & track_dese], col = "#6186DC", pch = 16, cex = .5)

mtext("latitude", side = 2, line = 9.1, las = 2, adj = 0)

addsubplotleg("(a)")


ylims <- range(unlist(xx))

# lon
plot(track_datenum, xx_m, type = 'n', axes = FALSE, xlab = '', ylab = '', ylim = ylims, xlim = xlims)
abline(h = c(-75, -74.8, -74.6), col = "grey85", lty = 2)

for(i in 1:100) {
  lines(track_datenum, xx[i, ], col = rgb(0, 0, 0, .02))
}

points(track_datenum[!track_filt & track_dese], xx_m[!track_filt & track_dese], pch = 16, cex = .5)
# axis(2, las = 1)
axis(2, las = 1, at = c(-75, -74.8, -74.6), labels = c("-75\u00B0 00\'", "-74\u00B0 48\'", "-74\u00B0 36\'"))
# legend("topleft", legend = "latitude", bty = 'n')
box(lwd = 1.5)
# rect(win[1], -10000, win[2], 10000, col = rgb(0, 0, 1, .25), border = NA)
abline(v = ceest, col = "#EF476F", lty = 2, lwd = 2)

points(track_datenum[track_filt & track_dese], xx_m[track_filt & track_dese], col = "#6186DC", pch = 16, cex = .5)
mtext("longitude", side = 2, line = 9.1, las = 2, adj = 0)

addsubplotleg("(b)")


# 
# # lon
# plot(track$datenum, track$X, type = 'n', axes = FALSE, xlab = '', ylab = '', ylim = c(-75.1, -74.5))
# abline(h = c(-75, -74.8, -74.6), col = "grey85", lty = 2)
# points(track$datenum, track$X, pch = 16, cex = .75)
# axis(2, las = 1, at = c(-75, -74.8, -74.6), labels = c("-75\u00B0 00\'", "-74\u00B0 48\'", "-74\u00B0 36\'"))
# # axis(2, las = 1)
# # legend("topleft", legend = "longitude", bty = 'n')
# box(lwd = 1.5)
# # rect(win[1], -10000, win[2], 10000, col = rgb(0, 0, 1, .25), border = NA)
# abline(v = ceest, col = "#EF476F", lty = 2, lwd = 2)
# points(track_filt$datenum, track_filt$X, col = "#6186DC", pch = 16, cex = .75)
# mtext("longitude", side = 2, line = 11.1, las = 2, adj = 0)


par(mar = c(3.1, 11.1, 0.1, 0.1))

xlims <- xlims/60/60 - min(ceest)/60/60

# big one
plot(stag$Date/60/60 - min(ceest)/60/60, -stag$Depth, type = 'n', axes = FALSE, xlab = '', ylab = '', ylim = c(-2050, 50), xlim = xlims)
abline(h = c(0, -400, -800, -1200, -1600, -2000), col = "grey85", lty = 2)
ustretchid <- unique(stag$stretchid)
for(i in 1:length(ustretchid)) {
  cur <- stag[stag$stretchid == ustretchid[i], ]
  lines(cur$Date/60/60 - min(ceest)/60/60, -cur$Depth)
}

# rect(win[1]/60/60 - min(ceest)/60/60, -10000, win[2]/60/60 - min(ceest)/60/60, 10000, col = rgb(0, 0, 1, .25), border = NA)
abline(v = ceest/60/60 - min(ceest)/60/60, col = "#EF476F", lty = 2, lwd = 2)
axis(2, las = 1, at = c(0, -400, -800, -1200, -1600, -2000), labels = c(0, 400, 800, 1200, 1600, 2000))
box(lwd = 1.5)
# x axis
axis(1)
# legend("topleft", legend = "stag", bty = 'n')
lines(stag_filt$Date/60/60 - min(ceest)/60/60, -stag_filt$Depth, type = 'l', col = "#6186DC")
mtext("depth (m)", side = 2, line = 9.1, las = 2, adj = 0)


addsubplotleg("(c)")

par(mar = c(0.1, 11.1, 0.1, 0.1))

# zoom
xlims <- c(win[1] - 0.25*60*60, win[2] + 0.25*60*60)
plot(stag_filt$Date, -stag_filt$Depth, type = 'n', axes = FALSE, xlab = '', ylab = '', xlim = xlims, ylim = c(-1600, 100))
abline(h = c(0, -400, -800, -1200), col = "grey85", lty = 2)
lines(stag$Date, -stag$Depth, lwd = 3)
box(col = "#6186DC", lwd = 1.5)
# axis(2, las = 1, col = "#6186DC")
axis(2, las = 1, col = "#6186DC", at = seq(0, -1400, by = -200), labels = abs(seq(0, -1400, by = -200)))
lines(dtag$t, -dtag$p, col = "#069D75", lwd = 2)
abline(v = ceest, col = "#EF476F", lty = 2, lwd = 2)
legend("bottomleft", legend = c("STag093", "DTag19_218", "CEE start"), col = c("black", "#069D75", "#EF476F"), lwd = c(3, 2, 2), lty = c(1, 1, 2), bty = 'n')
mtext("depth (m)", side = 2, line = 9.1, las = 2, adj = 0)


addsubplotleg("(d)")

par(mar = c(3.1, 11.1, 0.1, 0.1))



# calculate NSD from the cee
dis <- which(min(abs(ceest - track_datenum)) == (ceest - track_datenum))

nsd <- list()
for(i in 1:100) {
  nsd[[i]] <- latlond(yy[i, dis], xx[i, dis], yy[i, ], xx[i, ])^2
}

nsd_r <- do.call('rbind', nsd)
nsd_m <- colMeans(nsd_r)

ylims <- range(unlist(nsd_r[, track_datenum >= win[1] & track_datenum <= win[2]]))
xlims <- xlims/60/60 - min(ceest)/60/60
plot(track_datenum/60/60 - min(ceest)/60/60, nsd_m, type = 'n', axes = FALSE, xlab = '', ylab = '', xlim = xlims, ylim = ylims + c(-50, 50))
abline(h = seq(0, 400, by = 100), col = "grey85", lty = 2)

for(i in 1:100) {
  lines(track_datenum/60/60 - min(ceest)/60/60, nsd_r[i, ], col = rgb(0, 0, 0, .04))
}

# points(track_datenum[!track_dese & track_filt]/60/60 - min(ceest)/60/60, nsd_m[!track_dese & track_filt], cex = .75)
points(track_datenum[track_dese & track_filt]/60/60 - min(ceest)/60/60, nsd_m[track_dese & track_filt], pch = 16, cex = .75)


box(col = "#6186DC", lwd = 1.5)
axis(2, las = 1, col = "#6186DC")
abline(v = ceest/60/60 - min(ceest)/60/60, col = "#EF476F", lty = 2, lwd = 2)
# legend("topleft", legend = "net squared displacement", bty = 'n')
axis(1, col = "#6186DC")
mtext(expression(paste("NSD (", km^2, ")")), side = 2, line = 9.1, las = 1, adj = 0)


addsubplotleg("(e)")


par(mar = c(0.1, 11.1, 0.1, 0.1))




# dtag zoom

xlim_prime <- max(ceest - min(dtag$t), max(dtag$t) - ceest)
xlims <- c(ceest - xlim_prime, ceest + xlim_prime)/60/60 - min(ceest)/60/60

plot(dtag$t/60/60 - min(ceest)/60/60, -dtag$p, type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = c(-1500, 100), xlim = xlims)
abline(h = c(0, -400, -800, -1200), col = "grey85", lty = 2)
lines(dtag$t/60/60 - min(ceest)/60/60, -dtag$p, col = "black")
for(i in 1:nrow(dtag$dive_cues)) {
  lines(
    dtag$t[dtag$dive_cues$st[i]:dtag$dive_cues$en[i]]/60/60 - min(ceest)/60/60,
    -dtag$p[dtag$dive_cues$st[i]:dtag$dive_cues$en[i]],
    col = i + 1, lwd = 2
  )
}
box(col = "#069D75", lwd = 1.5)
axis(2, las = 1, col = "#069D75", at = seq(0, -1400, by = -200), labels = abs(seq(0, -1400, by = -200)))
# legend("topleft", legend = "dtag", bty = 'n')
abline(v = ceest/60/60 - min(ceest)/60/60, col = "#EF476F", lty = 2, lwd = 2)
mtext("depth (m)", side = 2, line = 9.1, las = 1, adj = 0)

text(rowMeans(cbind(dtag$t[dtag$dive_cues$st], dtag$t[dtag$dive_cues$en]))/60/60 - min(ceest)/60/60, -100, 1:nrow(dtag$dive_cues), col = 1:nrow(dtag$dive_cues) + 1)


addsubplotleg("(f)")

par(mar = c(4.1, 11.1, 0.1, 0.1))

# msa

msa <- 9.81 * abs(sqrt(rowSums(dtag$A^2)) - 1)


plot(dtag$t/60/60 - min(ceest)/60/60, msa, type = 'n', axes = FALSE, xlab = "", ylab = "", xlim = xlims)
# rect(dtag$t[dtag$dive_cues$st]/60/60 - min(ceest)/60/60, -20, dtag$t[dtag$dive_cues$en]/60/60 - min(ceest)/60/60, 20, col = 1:nrow(dtag$dive_cues) + 1, border = NA)
abline(h = c(0, 9.81/2, 9.81), col = "grey85", lty = 2)
lines(dtag$t/60/60 - min(ceest)/60/60, msa, lwd = 1)

for(i in 1:nrow(dtag$dive_cues)) {
  lines(
    dtag$t[dtag$dive_cues$st[i]:dtag$dive_cues$en[i]]/60/60 - min(ceest)/60/60,
    msa[dtag$dive_cues$st[i]:dtag$dive_cues$en[i]],
    col = i + 1, lwd = 1.5
  )
}

box(col = "#069D75", lwd = 1.5)
axis(2, las = 1, col = "#069D75", at = c(0, 4, 8))
axis(1, col = "#069D75")

# legend("topleft", legend = "pitch", bty = 'n')
abline(v = ceest/60/60 - min(ceest)/60/60, col = "#EF476F", lty = 2, lwd = 2)
mtext(expression(paste("MSA (", m %.% s^{-2}, ")")), side = 2, line = 9.1, las = 1, adj = 0)
mtext("hours since (to) sound exposure", side = 1, line = 3.1)
addsubplotleg("(g)")

dev.off()
