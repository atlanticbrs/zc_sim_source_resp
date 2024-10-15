register_target(
  # group the data by deployid
  
  tar_group_by(
    name = data_group,
    command = fit_data,
    deployid
  )
)