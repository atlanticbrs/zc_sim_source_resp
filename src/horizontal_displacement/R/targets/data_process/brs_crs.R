register_target(
  # Default Coordinate Reference System
  tar_target(
    name = brs_crs,
    command = {
      
      brs_crs <- paste("+proj=aea +lat_1=27.33333333333333",
                       "+lat_2=40.66666666666666 +lat_0=34 +lon_0=-78",
                       "+x_0=0 +y_0=0",
                       "+ellps=WGS84 +datum=WGS84 +units=m +no_defs")
      
      return(brs_crs)
    }
  )
)