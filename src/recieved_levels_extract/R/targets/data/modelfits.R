register_target(
  # fit the data with crawl
  
  tar_target(
    name = modelfits,
    command = {
      
      list(fit = crawlUtils::cu_crw_argos(data_group),
           deployid = unique(data_group$deployid))
      
    },
    pattern = map(data_group),
    iteration = "list"
  )
)