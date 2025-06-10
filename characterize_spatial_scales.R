#What are the relevant spatial scales associated with lion groups?

library(cocomo)

#PARAMETERS - MODIFY THESE

#directories
dir <- '~/Dropbox/lions_at_night/' #path to project directory
R <- 10 #radius (in meters) for spatial headings
speed_dt <- 10 #time difference to use for computing speeds and headings (in units of timesteps)
min_speed_to_compute_heading <- .5 #in m/s
start_times <- as.POSIXct(paste(seq.Date(date('2023-05-06'),date('2023-06-15'),by=1), '16:00:00'), tz = 'UTC') #timestamps of times to start each analysis period
end_times <- as.POSIXct(paste(seq.Date(date('2023-05-07'),date('2023-06-16'),by=1), '06:00:00'), tz = 'UTC') #timestamps of times to end each analysis period
dist_bins <- c(0, 10^seq(0.6,3.6,.4))

#------------------------------------------------
#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#FUNCTIONS

#LOAD DATA
gpsfile <- paste0(dir, 'data_namibia/processed/lion_omukutu_xy_level1.RData')
load(gpsfile) 

#basic info
n_inds <- nrow(xs)
n_times <- ncol(xs)

#GET PERIODS OF DATA TO USE
periods <- data.frame(start_time = start_times, end_time = end_times)
periods$t0_idx <- match(periods$start_time, timestamps)
periods$tf_idx <- match(periods$end_time, timestamps)

#relevant time indexes
relevant_t_idxs <- c()
for(i in 1:nrow(periods)){
  relevant_t_idxs <- c(relevant_t_idxs, periods$t0_idx[i]:periods$tf_idx[i])
}
idxs_to_remove <- setdiff(1:n_times, relevant_t_idxs)


#COMPUTE METRICS

#headings and speeds
heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
for(i in 1:n_inds){
  out <- cocomo::get_heading_and_speed_temporal(xs[i, ], ys[i,], t_window = speed_dt, forward = T)
  heads[i,] <- out$heads
  speeds[i,] <- out$speeds
  heads[i,which(speeds[i,] < min_speed_to_compute_heading)] <- NA
}

#remove heading when speed is too low

#dyadic distances
dyad_dists <- cocomo::get_group_dyadic_distances(xs, ys)

#dyadic distance changes, speed differences, and heading correlations between individuals
dyad_dist_changes <- head_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    head_corrs[i,j,] <- cos(heads[i,])*cos(heads[j,]) + sin(heads[i,])*sin(heads[j,])
    log_speed_diffs[i,j,] <- abs(log(speeds[i,]) - log(speeds[j,]))
    speed_diffs[i,j,] <- abs(speeds[i,] - speeds[j,])
    idxs_now <- seq(1,(n_times-speed_dt))
    idxs_fut <- seq(speed_dt + 1, n_times)
    dyad_dist_changes[i,j,idxs_now] <- dyad_dists[i,j,idxs_fut] - dyad_dists[i,j,idxs_now]
  }
}

#set non-relevant indexes to NA
dyad_dists[,,idxs_to_remove] <- NA
head_corrs[,,idxs_to_remove] <- NA
speed_diffs[,,idxs_to_remove] <- NA
log_speed_diffs[,,idxs_to_remove] <- NA
dyad_dist_changes[,,idxs_to_remove] <- NA

#calculate mean metrics during relevant periods only, as a function of distance apart
head_corrs_mean <- speed_diffs_mean <- log_speed_diffs_mean <- change_dyad_dist_mean <- rep(NA, length(dist_bins)-1)
for(i in 1:(length(dist_bins)-1)){
  idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
  head_corrs_mean[i] <- mean(head_corrs[idxs], na.rm=T)
  speed_diffs_mean[i] <- mean(speed_diffs[idxs], na.rm=T)
  log_speed_diffs_mean[i] <- mean(log_speed_diffs[idxs], na.rm=T)
  change_dyad_dist_mean[i] <- mean(dyad_dist_changes[idxs], na.rm=T)
}

#PLOTS

#Plot 1: Distribution of log(dyadic distances) between individuals
quartz()
histo <- hist(log(dyad_dists, 10), plot=T, breaks=50, xlab = 'Log (base 10) dyadic distance (m)',main='')

#Plot 2: Mean heading correlation as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],head_corrs_mean, xlab = 'Distance apart (m)', ylab = 'Mean heading correlation', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,1),log='x')

#Plot 3: Mean difference in speed as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],speed_diffs_mean, xlab = 'Distance apart (m)', ylab = 'Mean absolute speed difference (m/s)', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,max(speed_diffs_mean,na.rm=T)), log='x')

#Plot 4: Mean change in dyadic distance
quartz()
plot(dist_bins[2:length(dist_bins)],change_dyad_dist_mean, xlab = 'Distance apart (m)', ylab = 'Mean change in dyadic distance per timestep (m)', pch = 19, col = '#00000066', cex = 1.5, ylim = c(min(change_dyad_dist_mean,na.rm=T),max(change_dyad_dist_mean,na.rm=T)), log = 'x')
abline(h=0, lty = 2)

#TODO:
#show distribution (violion plot, median, IQR)
#look into autocorrelation, consider downsampling 
#model with autocorrelation model to get error bars 
