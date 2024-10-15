register_target(
  # Prepare the data to be fit by crawl
  tar_target(
    name = fit_data,
    command = {
      
      dat <- dplyr::mutate(locs_filter) %>%
        crawlUtils::cu_add_argos_cols()
      
      return(dat)
    }
  )
)