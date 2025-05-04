#Import high-res GPS data from csvs and convert to matrix form, preprocess to level 1

library(cocomo)
library(lubridate)
Sys.setenv(TZ='UTC')

#Directories
datadir <- '~/Dropbox/lions_at_night/data_namibia/rawdata/'
outdir <- '~/Dropbox/lions_at_night/data_namibia/processed/'
level0_outfile <- paste0(outdir, 'lion_omukutu_latlon_xy_level0.RData')
level1_outfile <- paste0(outdir, 'lion_omukutu_xy_level1.RData')
level1_outfile_latlon <- paste0(outdir, 'lion_omukutu_latlon_level1.RData')

#Read in data
alldat <- data.frame()
files <- list.files(datadir, full.names = T)
for(i in 1:length(files)){
  newdat <- read.csv(files[i],stringsAsFactors = F, header=T)
  alldat <- rbind(alldat, newdat)
}

#rename to correct column names
#timestamps are UTC
colnames(alldat) <- c('individual.local.identifier','Timestamp','location.lat','location.long','height')

#convert timestamp to POSIXct
alldat$timestamp <- lubridate::as_datetime(alldat$Timestamp, tz = 'UTC', format = '%d/%m/%Y %H:%M:%OS')

#make data chunks data frame
data_chunks <- data.frame(start = min(alldat$timestamp), end = max(alldat$timestamp))

#convert to matrix form and save output file
cocomo::reformat_movebank_to_matrix(movebank_data = alldat, 
                                    output_file_path = level0_outfile, 
                                    data_chunks = data_chunks, 
                                    utm_zone = 33, 
                                    hemisphere = 'south', 
                                    output_utm = T, 
                                    output_latlon = T)

#load the level 0 file we just created
load(level0_outfile)

#preprocess from level 0 to level 1
out <- preprocess_gps_level0_to_level1(input_file_path = level0_outfile, 
                                output_file_path = level1_outfile, 
                                xs = xs, 
                                ys = ys, 
                                timestamps = timestamps, 
                                ids = ids, 
                                remove_unrealistic_speeds = T, 
                                remove_isolated_points = T, 
                                remove_unrealistic_locations = T, 
                                interpolate_small_gaps = T, 
                                interpolate_stationary_periods = T)

#also get lat/lon data at level 1
#get lon/lat data (level 1)
print('converting back to lat/lon to save latlon level 1 file')
lons <- lats <- matrix(nrow = nrow(xs), ncol = ncol(out$xs))
for(i in 1:nrow(xs)){
  lonsLats <- cocomo::utm_to_latlon(cbind(out$xs[i,],out$ys[i,]),utm_zone = 33, hemisphere = 'south')
  lons[i,] <- lonsLats[,1]
  lats[i,] <- lonsLats[,2]
}

save(file = level1_outfile_latlon, list = c('lats','lons','timestamps','ids'))
