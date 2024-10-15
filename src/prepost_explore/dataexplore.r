###
# quick data exploration of the pre/post results
#

# required libraries
library(corrplot)

# load data
dat <- readRDS("full_prepost_results.rds")
meta <- read.table("../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# rectify differences in deployid conventions
uids <- unique(meta$deployid)
strippeduids <- sub("_DUML", "", uids)

desetofix <- dat$tag %in% strippeduids[grep("_DUML", uids)]
dat$tag[desetofix] <- paste0(dat$tag[desetofix], "_DUML")

# filter exposure events to just ones with control, scaled source with series data
meta <- meta[meta$cee_type %in% c("control", "simulated mfas"), ]
meta <- meta[meta$ser, ]

# reshape into wide
dat$lab <- paste(dat$tag, dat$cee_idx)
dat$statlab <- paste(dat$stat, paste0("L=", dat$lag_len), paste0("W=", dat$window_len))
dat2 <- dat[, c('lab', 'statlab', 'p.left')]
dat3 <- reshape(dat2, idvar = 'lab', timevar = 'statlab', direction = 'wide')
rownames(dat3) <- dat3$lab
dat4 <- dat3[, -1]

# now do the right p-value
dat2 <- dat[, c('lab', 'statlab', 'p.right')]
dat3 <- reshape(dat2, idvar = 'lab', timevar = 'statlab', direction = 'wide')
rownames(dat3) <- dat3$lab
dat4right <- dat3[, -1]

# try to match back up
tagid_ceeidx <- strsplit(rownames(dat4), ' ')
tagid <- sapply(tagid_ceeidx, '[[', 1)
ceeidx <- sapply(tagid_ceeidx, '[[', 2)
ceeid <- vector()

for(i in 1:length(unique(tagid))) {
  curid <- which(tagid == unique(tagid)[i])
  curmeta <- meta[meta$deployid == tagid[curid[1]], ]
  curmeta <- curmeta[curmeta$cee_type %in% c("simulated mfas", "control"), ]
  curcee <- sort(curmeta$cee_id)[as.numeric(ceeidx[curid])]
  ceeid <- c(ceeid, curcee)
}

# get rid of anything without a cee
dat4 <- dat4[!is.na(ceeid), ]
dat4right <- dat4right[!is.na(ceeid), ]

meta$key <- paste(meta$deployid, meta$cee_id)

# remake row labels
dat4_deployid <- sapply(strsplit(rownames(dat4), " "), '[[', 1)
dat4_ceeid <- ceeid[!is.na(ceeid)]
dat4_rlbin <- meta$rl_bin[match(paste(dat4_deployid, dat4_ceeid), meta$key)]
dat4_rlbin <- factor(dat4_rlbin, levels = c("control", "<100", "100-120", ">120"))

rownames(dat4) <- paste(sub("_DUML", "", sub("Zc", "", dat4_deployid)), dat4_ceeid, dat4_rlbin)
rownames(dat4right) <- paste(sub("_DUML", "", sub("Zc", "", dat4_deployid)), dat4_ceeid, dat4_rlbin)
colnames(dat4) <- sapply(strsplit(colnames(dat4), "\\."), '[[', 3)
colnames(dat4right) <- sapply(strsplit(colnames(dat4right), "\\."), '[[', 3)

# correlations
# datcor <- cor(dat4)
# corrplot::corrplot(datcor, method = 'color', order = "hclust", tl.cex = 0.5, main = 'left p')
# datcor <- cor(dat4right)
# corrplot::corrplot(datcor, method = 'color', order = "hclust", tl.cex = 0.5, main = 'right p')

# corrplot::corrplot(as.matrix(dat4), is.corr = FALSE, method = 'color', col.lim = c(0, 1), col = COL2('RdBu', 200), tl.cex = 0.5, main = 'left p')
# corrplot::corrplot(as.matrix(dat4right), is.corr = FALSE, method = 'color', col.lim = c(0, 1), col = COL2('RdBu', 200), tl.cex = 0.5, main = 'right p')

# reorder by most common
sigmat <- as.matrix(dat4)
sigmat <- sigmat > 0.95 | sigmat < 0.05

oo <- order(colSums(sigmat), decreasing = TRUE)
mdato <- as.matrix(dat4)[, oo]

oo2 <- order(rowSums(sigmat), decreasing = TRUE)
mdato2 <- mdato[oo2, ]

oo3 <- order(dat4_rlbin[oo2], decreasing = TRUE)
mdato3 <- mdato2[oo3, ]

# corrplot::corrplot(mdato3, is.corr = FALSE, method = 'color', col.lim = c(0, 1), col = COL2('RdBu', 200), tl.cex = 0.5, main = 'all metric sorted')



### collapse the metrics
metrics_short <- sapply(strsplit(colnames(dat4), ' '), '[[', 1)
metrics_short <- sub("timing_", "", metrics_short)

umetrics_short <- unique(metrics_short)
nmetrics_short <- length(umetrics_short)

# collapse metrics for left
collapsed_metrics_left <- list()
for(i in 1:nmetrics_short) {
  cur <- metrics_short == umetrics_short[i]
  
  met_cur <- do.call('rbind', apply(dat4[, cur], 1, function(x) {
    data.frame(
      siglo = length(which(x < 0.05)) / length(which(cur))#,
      # sighi = length(which(x > 0.95)) / length(which(cur))
    )
  }))
  
  names(met_cur) <- paste(umetrics_short[i], c("low"))
  collapsed_metrics_left[[i]] <- met_cur
}

# collapse metrics for right
collapsed_metrics_right <- list()
for(i in 1:nmetrics_short) {
  cur <- metrics_short == umetrics_short[i]
  
  met_cur <- do.call('rbind', apply(dat4right[, cur], 1, function(x) {
    data.frame(
      siglo = length(which(x < 0.05)) / length(which(cur))#,
      # sighi = length(which(x > 0.95)) / length(which(cur))
    )
  }))
  
  names(met_cur) <- paste(umetrics_short[i], c("high"))
  collapsed_metrics_right[[i]] <- met_cur
}

collapsed_metrics_left <- do.call(cbind, collapsed_metrics_left)
collapsed_metrics_right <- do.call(cbind, collapsed_metrics_right)
collapsed_metrics <- cbind(collapsed_metrics_left, collapsed_metrics_right)

oo_col <- order(colSums(collapsed_metrics), decreasing = TRUE)
oo_row <- order(rowSums(collapsed_metrics), decreasing = TRUE)
oo_bin <- order(dat4_rlbin[oo_row], decreasing = TRUE)

# make the mar large enough for the angled labels

# pdf("test1.pdf", width = 10, height = 10)
# corrplot::corrplot(as.matrix(collapsed_metrics[oo_row[oo_bin], oo_col]), is.corr = FALSE, method = 'color', col.lim = c(0, 0.89), tl.cex = 1, tl.col = "black", tl.srt = 60, na.label = "/", cl.pos = 'b')
# dev.off()



### make one where high and low are separated?
collapsed_metrics_lo <- collapsed_metrics[, grep('low', colnames(collapsed_metrics))]
collapsed_metrics_hi <- collapsed_metrics[, grep('high', colnames(collapsed_metrics))]

oo_col_hi <- order(colSums(collapsed_metrics_hi), decreasing = TRUE)
oo_col_lo <- order(colSums(collapsed_metrics_lo), decreasing = TRUE)

# par(mfrow = c(1, 2))
# corrplot::corrplot(as.matrix(collapsed_metrics_hi[oo_row[oo_bin], oo_col_hi]), is.corr = FALSE, method = 'color', tl.cex = 1, tl.col = "black", tl.srt = 60, na.label = "/", cl.pos = 'b')
# corrplot::corrplot(as.matrix(collapsed_metrics_lo[oo_row[oo_bin], oo_col_lo]), is.corr = FALSE, method = 'color', tl.cex = 1, tl.col = "black", tl.srt = 60, na.label = "/", cl.pos = 'b')



### make it all one matrix with 0 colum in the middle
blocked_collapsed_metrics <- cbind(collapsed_metrics_hi, cbind(rep(NA, nrow(collapsed_metrics))), collapsed_metrics_lo)
names(blocked_collapsed_metrics)[11] <- " "

# fix other labels
rn <- rownames(blocked_collapsed_metrics)
rn <- sub('120', '120 dB SPL', rn)
rn <- sub('<100', '<100 db SPL', rn)
rownames(blocked_collapsed_metrics) <- sub('Tag', 'STag', rn)
cn <- colnames(blocked_collapsed_metrics)
cn <- sub(' low', '', cn)
cn <- sub(' high', '', cn)
cn <- sub('avg', 'mean', cn)
cn <- gsub('_', ' ', cn)
colnames(blocked_collapsed_metrics) <- cn

blocked_collapsed_metrics <- blocked_collapsed_metrics[oo_row[oo_bin], ]

rownames(blocked_collapsed_metrics)[grep(">120", rownames(blocked_collapsed_metrics))[1]] <- paste(   ">120 dB SPL     ", rownames(blocked_collapsed_metrics)[grep(">120", rownames(blocked_collapsed_metrics))[1]])
rownames(blocked_collapsed_metrics)[grep("100-120", rownames(blocked_collapsed_metrics))[1]] <- paste("100-120 dB SPL  ", rownames(blocked_collapsed_metrics)[grep("100-120", rownames(blocked_collapsed_metrics))[1]])
rownames(blocked_collapsed_metrics)[grep("<100", rownames(blocked_collapsed_metrics))[1]] <- paste(   "<100 dB SPL     ", rownames(blocked_collapsed_metrics)[grep("<100", rownames(blocked_collapsed_metrics))[1]])
rownames(blocked_collapsed_metrics)[grep("control", rownames(blocked_collapsed_metrics))[1]] <- paste("control         ", rownames(blocked_collapsed_metrics)[grep("control", rownames(blocked_collapsed_metrics))[1]])

rownames(blocked_collapsed_metrics) <- gsub(">120 dB SPL$", "", rownames(blocked_collapsed_metrics))
rownames(blocked_collapsed_metrics) <- gsub("100-120 dB SPL$", "", rownames(blocked_collapsed_metrics))
rownames(blocked_collapsed_metrics) <- gsub("<100 db SPL$", "", rownames(blocked_collapsed_metrics))
rownames(blocked_collapsed_metrics) <- gsub("control$", "", rownames(blocked_collapsed_metrics))



### make plot
pdf("prepostexplore.pdf", width = 8, height = 7)
par(mar = c(0, 0, 5.1, 0), oma = c(0, 0, 0, 0), family = "serif")
corrplot::corrplot(as.matrix(blocked_collapsed_metrics), is.corr = FALSE, method = 'color', tl.cex = 1, tl.col = "black", tl.srt = 60, na.label = " ", cl.pos = 'r', col = corrplot::COL1("Blues"))

xseq <- seq(0.5, 0.5 + 10, by = 1)
yseq <- seq(0.5, 0.5 + nrow(blocked_collapsed_metrics), by = 1)
segments(xseq, yseq[1], xseq, yseq[length(yseq)], col = "grey85", lwd = 2)
segments(xseq[1], yseq, xseq[length(xseq)], yseq, col = "grey85", lwd = 2)

xseq <- seq(0.5 + 11, 0.5 + 21, by = 1)
yseq <- seq(0.5, 0.5 + nrow(blocked_collapsed_metrics), by = 1)
segments(xseq, yseq[1], xseq, yseq[length(yseq)], col = "grey85", lwd = 2)
segments(xseq[1], yseq, xseq[length(xseq)], yseq, col = "grey85", lwd = 2)

par(xpd = TRUE)
mtext('right tail', 3, line = 4.1, at = 7)
mtext('left tail', 3, line = 4.1, at = 18)
par(xpd = FALSE)

dev.off()
