###
# output univariate plots for supplement
# ~ wrc

# load data
plot_inputs <- readRDS("../02_intermediate_outputs/plot_inputs.rds")

# make a key for the labels
whatkey <- data.frame(key = unique(sapply(plot_inputs, '[[', 'what')), lab = NA)
whatkey[, 2] <- c(
  "IDDI (minutes)",
  "prior IDDI (minutes)",
  "dive duration (minutes)",
  "dive max. depth (meters)",
  "dive ascent duration (minutes)",
  "no. bounce dives",
  "mean bounce dive duration (minutes)",
  "mean bounce dive max. depth (meters)",
  "IDDI max. depth (meters)"
)

make_metric_histo <- function(cur, ...) {
  obs <- cur$obs
  bas <- cur$bas
  
  # make xlim to include the exp dive
  xlims <- range(c(obs, bas), na.rm = TRUE) * c(0.8, 1.2) # 20% slop
  
  xlabel <- whatkey$lab[which(cur$what == whatkey$key)]
  
  # change bounce label
  xlabel <- gsub ("bounce", "shallower dive", xlabel)
  
  hh <- hist(bas,
             # xlim = xlims,
             main = "",
             xlab = xlabel,
             ylab = "frequency",
             yaxt = 'n',
             xlim = xlims,
             ...
  )
  
  axis(2, las = 1)
  
  # make obs arrow
  arrowheight <- max(hh$counts) * 0.40
  arrows(obs, arrowheight, obs, 0, length = 0.15, col = "purple", lwd = 2)
  
  # make labels
  bas_nona <- bas[!is.na(bas)]
  perc = round(length(which(obs > bas_nona)) / length(bas_nona), 3)
  legend("topright", 
    c(cur$info$label_selected, 
      paste("prog =", cur$info$settings), 
      paste("quant =", round(perc, 2))
    ),
    bty = 'n'
  )
  
  box()
}


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

### meta

# load rls and make a key
rls <- read.table("../../../01_shared_data_products/event_manifest/exposure_events.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# filter
# rls$rl_bin_spl[rls$median_range_km < MAX_MODELLING_DISTANCE_RL_KM & rls$ceetype == 'control'] <- 'control'
rls <- rls[!is.na(rls$rl_bin_spl), ]
rlkey_rls <- paste(rls$deployid, rls$cee_id)


### inputs plots
plot_inputs <- readRDS("../02_intermediate_outputs/plot_inputs.rds")
plot_type <- sapply(plot_inputs, '[[', 'what')
deployid <- sapply(plot_inputs, function(x) x$info$deployid)

udeployid <- unique(deployid)
count <- 1
caption <- vector()
for(i in 1:length(udeployid)) {
  dose <- which(deployid == udeployid[i])
  cees <- sapply(plot_inputs[dose], function(x) x$info$ceeid)
  ucees <- unique(cees)

  for(q in 1:length(ucees)) {
  dese <- dose[which(cees == ucees[q])]
  
  ht <- 7 * 320
  wd <- 7 * 320
  mfrows <- c(3, 3)
  if(length(dese) == 2) {
    ht = 10 * 320
    wd = 5 * 320
    mfrows <- c(2, 1)
  }
  
  deployid_label <- sub("Zc", "S", plot_inputs[[dese[1]]]$info$deployid)
  deployid_label <- sub("_DUML", "", deployid_label)
    
  key <- paste(plot_inputs[[dese[1]]]$info$deployid, plot_inputs[[dese[1]]]$info$ceeid)
  rlbin <- rls$rl_bin_spl[which(rlkey_rls == key)]
if(length(rlbin) > 0) {
  png(paste0(deployid_label, "-", plot_inputs[[dese[1]]]$info$ceeid, ".png"), height = ht, width = wd, res = 320)
  par(mar = c(4.1, 3.1, 3.1, 1.1), mfrow = mfrows, family = 'serif')
  
  caption[i] <- paste0("Figure S", i+37, ". ", deployid_label, " dive metric distributions for exposure trial ", plot_inputs[[dese[1]]]$info$ceeid)
  count <- count + 1
  
  if(rlbin != "control") rlbin <- paste(rlbin, "dB SPL")
  title_label <- paste(deployid_label, plot_inputs[[dese[1]]]$info$ceeid, rlbin)
  
  for(p in 1:length(dese)) {
    cur <- plot_inputs[[dese[p]]]
    make_metric_histo(cur)
    if(mfrows[1] == 2 & p == 1) title(title_label)
    if(mfrows[1] == 3 & p == 2) title(title_label)
  }
  dev.off()
}
}
}

write.table(caption, file = "captions.csv", row.names = FALSE, sep = ',')
