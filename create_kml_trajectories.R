#Create KMLs to animate GPS data in xs, ys, timestamps, ids format

#PARAMETERS
icon_file <- 'red_dot.png'
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

cols <- rep('white',n_inds)

#get dates
dates <- unique(date(timestamps))

#make a kml for each day
for(day in 1:length(dates)){
  
  #get beginning and end of day
  t0 <- (day-1)*60*60*24 + start_hr_utc*60*60 + 1
  tf <- t0 + 24*60*60 - 1
  
  filename <- paste0('~/Dropbox/lions_at_night/viz/kmls_omukutu/omukutu_',dates[day],'.kml')
  
  #get data
  lons_curr <- lons[,seq(t0,tf,step)]
  lats_curr <- lats[,seq(t0,tf,step)]
  timestamps_curr <- timestamps[seq(t0,tf,step)]
  
  # START WRITING
  sink(filename)
  
  # start output
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  cat("<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n")
  cat("<Document>\n")
  
  #locations - perhaps we want to uncomment and add some in
  # for(i in 1:nrow(den_latlon)){
  #   cat("<Placemark>\n")
  #   cat(paste("<name>",den_names[i],"</name>\n", sep = ' '))
  #   cat("<Point>\n")
  #   cat(paste("<coordinates>",den_latlon[i,1],den_latlon[i,2],'</coordinates>',sep=' '))
  #   cat("</Point>\n")
  #   cat(paste("<Icon><href>",den_icon,"</href></Icon>"))
  #   cat("</Placemark>\n")
  # }
  
  #tracks of ids
  for (i in 1:n_inds) {
    cat(paste("<Style id=\"track-",ids$code[i],"\"><IconStyle><scale>0.5</scale><Icon><href>",ids$icon[i],"</href></Icon></IconStyle><LineStyle><color>",as.character(cols[i]),"</color><colorMode>normal</colorMode></LineStyle></Style>\n",sep=""))
  }
  
  #time strings
  timestamps_curr <- as.character(format(timestamps_curr, '%y-%m-%d %H:%M:%S'))
  timestamps_curr <- gsub(' ','T',timestamps_curr)
  timestamps_curr <- paste0(timestamps_curr,'.000Z')
  
  #locations 
  for (i in 1:n_inds) {
    cat("<Folder>\n")
    cat(paste("<name>",ids$code[i],"</name>\n",sep=""))
    cat("<Placemark>")
    cat(paste("<styleUrl>#track-",ids$code[i],"</styleUrl>\n",sep=""))
    cat("<visibility>0</visibility>\n")
    cat(paste("<name>",ids$code[i],"</name>\n",sep=''))
    cat("<gx:Track>\n")
    cat("<gx:altitudeMode>relativeToGround</gx:altitudeMode>\n")
    
    # FOR EACH TIME
    for (tt in 1:length(timestamps_curr)) {
      if(!is.na(lons_curr[i,tt])){
        cat(sprintf("<when>%s</when>\n",timestamps_curr[tt]),sep="")
      } 
    }
    # FOR EACH TIME
    for (tt in 1:length(timestamps_curr)) {
      if(!is.na(lons_curr[i,tt])){
        cat(sprintf("<gx:coord>%s</gx:coord>\n",paste(lons_curr[i,tt],lats_curr[i,tt])),sep="")
      }
    }
    cat("</gx:Track>\n")
    cat("</Placemark>\n")
    cat("</Folder>\n")
  }
  
  cat("</Document>\n")
  cat("</kml>\n")
  sink()
}