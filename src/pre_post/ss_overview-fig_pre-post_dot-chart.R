# Make Panel Plot
# goal is to make a plot summarizing the pre-post results
library(ggplot2)
library(tidyverse)
library(tidytext)
library(extrafont)
# Pre-post results
df_all <- readr::read_csv('output/scaled_source_24hr/2024-03-04_pre-post_multiple-exposure.csv') %>% 
  dplyr::filter(tag != "ZcTag109")

# Metadata for a) labeling, and b) making sure we don't get multiple vlaues for ZcTag093
tags_mult_bline_metadata <- read_csv(file = 'data/sattag/tag_info_mult-segments_baseline_metadata.csv') %>% 
  dplyr::select(deployid, cee_id) %>% 
  distinct()

# RLs
rl_vals <- readr::read_csv('/Users/rob/Documents/_research/rrr/scaled_source-analysis-register-targets/results/output_data/rls_combined.csv') 
idx <- which(rl_vals$deployid == 'ZcTag093' & rl_vals$cee_id == '19_04')
rl_vals <- rl_vals[-idx, ]
rl_check <- inner_join(rl_vals, tags_mult_bline_metadata, by = c('deployid', 'cee_id')) %>% 
  mutate(new_id = paste(str_replace(deployid, 'Zc', 'S'), cee_id))

rl_vals$tag <- rl_vals$deployid
idx <- which(rl_vals$tag %in% c('ZcTag073', 'ZcTag080')) 
rl_vals$rl_bin[idx] <- 'Control'
rl_vals <- rl_vals %>% 
  mutate(new_id = paste(str_replace(deployid, 'Zc', 'S'), cee_id))

# Join in the rl_bin
df_joined <- df_all %>%
  left_join(rl_vals %>% select(tag, rl_bin, new_id), by = "tag") 

idx1 <- which(df_joined$tag == 'ZcTag069_cee-1802')
idx2 <- which(df_joined$tag == 'ZcTag069_cee-1803')
df_joined$tag[idx1] <- "STag069 18_02"
df_joined$tag[idx2] <- "STag069 18_03"
df_joined$rl_bin[idx1] <- "100-120"
df_joined$rl_bin[idx2] <- ">120"
df_joined$new_id[idx1] <- "STag069 18_02"
df_joined$new_id[idx2] <- "STag069 18_03"


df_joined <- df_joined %>% 
  mutate(rl_bin = replace_na(rl_bin, "Control"))

# Update last Missing IDs
miss_ids <- df_joined$tag[is.na(df_joined$new_id)]
new_cees <- tags_mult_bline_metadata %>% 
  filter(deployid %in% miss_ids)
new_ids <- paste(str_replace(miss_ids, 'Zc', 'S'), c('18_07', '22_02', '22_02'))
df_joined$new_id[is.na(df_joined$new_id)] <- new_ids

# change levels for plotting
df_joined <- df_joined %>% 
  mutate(rl_bin = recode(rl_bin, ">120" = ">120 dB SPL", 
                         "100-120" = "100-120 dB SPL", 
                         "<100" = "<100 dB SPL",
                         "Control" = "Control"))


# make a plot
# chance number of significances
nchance <- 1
# ggplot(df_all, aes(x = excess_signif, y = tag))+
p_strip <- ggplot(df_joined, aes(x = excess_signif, y = reorder_within(new_id, excess_signif, rl_bin)) )+
  scale_y_reordered() +
  # geom_vline(xintercept = nchance)+
  # geom_segment(aes(yend = tag), xend = 0, colour = 'grey50')+
  geom_point(size = 3)+
  # facet_grid(bin ~., scales = 'free')+
  # facet_grid(vars(factor(rl_bin, levels = c(">120",  "100-120", "<100",  "Control"))),
  facet_grid(vars(factor(rl_bin, levels = c(">120 dB SPL",  "100-120 dB SPL", "<100 dB SPL",  "Control"))),
             scales = 'free')+
  labs(y = 'Deploy ID', x = 'Odds Ratio of Excess Changes')+
  theme_bw(base_size = 12)+
  theme(text = element_text(family = "Times"))
             


