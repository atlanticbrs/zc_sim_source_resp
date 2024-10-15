file_path <- list.files('output/scaled_source_24hr', 
                        pattern="*.rds", full.names=TRUE)

library(tidyverse)
library(gt)

## counting how many files in your directory
file_number <- 1 

for (file in file_path) {
  name <- paste0("df_", file_number) ## getting a new name ready
  x <- readRDS(file) ## reading in the file
  assign(name, x) ## assigning the data the name we created two lines up
  file_number <- file_number + 1 
}

full_data <- mget(ls(pattern = "^df_")) %>%
  reduce(bind_rows)

rm(list = ls(pattern = glob2rx("df_*")))
saveRDS(full_data, 'output/scaled_source_24hr/full_prepost_results.rds')

df <- full_data %>% 
  group_by(tag) %>% 
  summarise(signif_level = 0.05,
            total_signif = sum(p.left < signif_level | p.right < signif_level), 
            prop_signif = total_signif / n(),
            total_signif_chance = n() * signif_level,
            excess_signif = total_signif / total_signif_chance) %>% 
  arrange(excess_signif)
readr::write_csv(df, 'output/scaled_source_24hr/2023-08-25_pre-post_single-exposure.csv')

gt(df) %>% 
  fmt_number(columns = c(excess_signif, prop_signif), 
             decimals = 2) 

df_control <- df %>% 
  filter(tag %in% c("ZcTag073", "ZcTag075", "ZcTag080", "ZcTag131", "ZcTag132"))%>% 
  arrange(tag)

gt(df_control) %>% 
  fmt_number(columns = c(excess_signif, prop_signif), 
             decimals = 2) 

# To do this next chunk, run the prepost_function_multiple-baseline.R script
# ZcTag069
df_mult <- full_data %>% 
  filter(!(tag == 'ZcTag131' & cee_idx == 1),
         !(tag == 'ZcTag132' & cee_idx == 1)) %>% 
  group_by(tag, cee_idx) %>% 
  summarise(signif_level = 0.05,
            total_signif = sum(p.left < signif_level | p.right < signif_level), 
            prop_signif = total_signif / n(),
            total_signif_chance = n() * signif_level,
            excess_signif = total_signif / total_signif_chance) %>% 
  arrange(excess_signif)
readr::write_csv(df_mult, 'output/scaled_source_24hr/2023-08-25_pre-post_multiple-exposure.csv')

gt(df_mult) %>% 
  fmt_number(columns = c(excess_signif, prop_signif), 
             decimals = 2)

# prep for binding
df_single <- df_mult %>% 
  filter(tag != "ZcTag069")

df_tag069 <- df_mult %>% 
  filter(tag == "ZcTag069") %>% 
  mutate(tag = c('ZcTag069_cee-1803', 
                 'ZcTag069_cee-1802')) 

df_all <- bind_rows(df_single, df_tag069)%>% 
  select(tag, signif_level,
         total_signif, prop_signif, 
         total_signif_chance, excess_signif) 
readr::write_csv(df_all, 'output/scaled_source_24hr/2024-03-04_pre-post_multiple-exposure.csv')

