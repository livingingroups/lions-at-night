#Create KMLs to animate GPS data in xs, ys, timestamps, ids format
library(cocomo)
library(lubridate)


#PARAMETERS
step <- 2 #resolution
start_hr_utc <- 12 #time to start animations for each day

#MAIN
#Load data
load('~/Dropbox/lions_at_night/data_namibia/processed/lion_omukutu_latlon_level1.RData')
n_inds <- nrow(ids)

#icon colors - blue = subadult male, red = adult female, orange = subadult females
ids$icon <- 'orange_dot.png'
ids$icon[which(ids$code=='OPL26')] <- 'blue_dot.png'
ids$icon[which(ids$code=='OPL29')] <- 'red_dot.png'

cols <- rep('ff377ef7',n_inds)
cols[which(ids$code=='OPL26')] <- 'ffed8031'
cols[which(ids$code=='OPL29')] <- 'ff0000ff'

#get dates
dates <- unique(lubridate::date(timestamps))

#make a kml for each day
for(day in 1:(length(dates)-1)){
  
  #get beginning and end of day
  t0 <- (day-1)*60*60*24 + start_hr_utc*60*60 + 1
  tf <- t0 + 24*60*60 - 1
  
  filename <- paste0('~/Dropbox/lions_at_night/viz/kmls_omukutu/omukutu_',dates[day],'.kml')
  
  cocomo::create_trajectories_kml(lons = lons, lats = lats, timestamps = timestamps, 
                                  id_codes = ids$code, t0 = t0, tf = tf, output_file_path = filename, 
                                  step = step, cols = cols, icons = ids$icon)
  
}