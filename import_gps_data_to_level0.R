#Import high-res GPS data from csvs and convert to matrix form, preprocess to level 1

library(cocomo)
library(lubridate)
Sys.setenv(TZ='UTC')

#Directories
datadir <- '~/Dropbox/lions_at_night/data_namibia/rawdata/'
outdir <- '~/Dropbox/lions_at_night/data_namibia/processed/'
outfile <- paste0(outdir, '/lion_omukutu_latlon_xy_level0.RData')

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
cocomo::reformat_movebank_to_matrix(movebank_data = alldat, output_file_path = outfile, data_chunks = data_chunks, utm_zone = 33, hemisphere = 'south', output_utm = T, output_latlon = T)
