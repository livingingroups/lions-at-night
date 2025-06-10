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
dist_bins <- c(0, 10^seq(0.6,4,.4))

#------------------------------------------------
#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#FUNCTIONS

#helper function to get statistics of the distribution of values
get_distrib_statistics <- function(vals, bins){
  out <- list()
  out$bins <- bins
  out$mean <- out$median <- out$q75 <- out$q25 <- out$q025 <- out$q975 <- rep(NA, length(bins)-1)
  #loop over all bins
  for(i in 1:(length(bins)-1)){
    
    #get indexes to data in that bin
    idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
    
    #get relevant stats
    out$mean[i] <- mean(vals[idxs], na.rm=T)
    out$median[i] <- median(vals[idxs], na.rm=T)
    out$q75[i] <- quantile(vals[idxs], 0.75, na.rm=T)
    out$q25[i] <- quantile(vals[idxs], 0.25, na.rm=T)
    out$q025[i] <- quantile(vals[idxs], 0.025, na.rm=T)
    out$q975[i] <- quantile(vals[idxs], 0.975, na.rm=T)
  }
  
  return(out)
}

#helper function to make plot
make_plot <- function(dat, plotpath, ylab, ylim = NULL, abline_y = NULL){
  if(is.null(ylim)){
    ylim <- c(min(dat$q25,na.rm=T), max(dat$q75, na.rm=T))
  }
  
  png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
  mids <- (dat$bins[2:length(dat$bins)] + dat$bins[1:(length(dat$bins)-1)]) / 2
  plot(mids, dat$mean, xlab = 'Distance apart (m)', ylab = ylab, pch = 19, col = '#00000066', cex = 1.5, ylim = ylim ,log='x')
  if(!is.null(abline_y)){
    abline(h=abline_y, lty = 2)
  }
  arrows(mids, dat$q25, mids, dat$q75, lwd = 2, code = 3, length = 0.1, angle = 90)
  points(mids, dat$median, pch = 19, col = 'black',cex = 1.5)
  points(mids, dat$mean, pch = 19, col = 'gray', cex = 1.5)
  dev.off()
}

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
dyad_dist_changes <- head_corrs <- speed_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
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

#calculate mean, median, iqr of metrics during relevant periods only, as a function of distance apart
metrics <- list()
metrics$dist_bins <- dist_bins

#get statistics of the distribution for each metric within each bin
metrics$head_corr <- get_distrib_statistics(head_corrs, dist_bins)
metrics$speed_diff <- get_distrib_statistics(speed_diffs, dist_bins)
metrics$log_speed_diff <- get_distrib_statistics(log_speed_diffs, dist_bins)
metrics$dyad_dist_change <- get_distrib_statistics(dyad_dist_changes, dist_bins)
metrics$approach <- get_distrib_statistics(dyad_dist_changes[which(dyad_dist_changes!=0)] < 0, dist_bins)
metrics$ang_between_heads <- get_distrib_statistics(acos(head_corrs)*180/pi, dist_bins)

#PLOTS

#Plot 1: Distribution of log(dyadic distances) between individuals
png(filename = paste0(dir, 'plots/spatial_scales/hist_log_dyad_dist.png'), width = 8, height = 6, units = 'in', res = 300)
histo <- hist(log(dyad_dists, 10), plot=T, breaks=80, xlab = 'Dyadic distance (m) - log bins',main='', xlim = c(0,4), freq = F, xaxt ='n')
axis(1, at = seq(0,4,1), labels = 10^seq(0,4,1))
dev.off()

#Plot 2: Angle between headings 
plotpath <- paste0(dir, 'plots/spatial_scales/ang_between_heads.png')
ylab <- "Angle between headings (degrees)"
ylim <- c(0,180)
dat <- metrics$ang_between_heads
make_plot(dat, plotpath, ylab, ylim, abline_y = 90)

#Plot 3: Heading correlation
plotpath <- paste0(dir, 'plots/spatial_scales/heading_corr.png')
ylab <- "Heading correlation"
ylim <- c(-1,1)
dat <- metrics$head_corr
make_plot(dat, plotpath, ylab, ylim, abline_y = 0)

#Plot 4: Mean difference in speed as a function of distance 
plotpath <- paste0(dir, 'plots/spatial_scales/abs_speed_diff.png')
ylab = 'Absolute speed difference (m/s)'
dat <- metrics$speed_diff
make_plot(dat, plotpath, ylab)

#Plot 5: Change in dyadic distance
plotpath <- paste0(dir, 'plots/spatial_scales/dyad_dist_change.png')
ylab = 'Change in dyadic distance (m)'
dat <- metrics$dyad_dist_change
make_plot(dat, plotpath, ylab, abline_y = 0)

#Plot 6: Probability of approaching - needs a different approach to error bars etc
#plotpath <- paste0(dir, 'plots/spatial_scales/p_approach.png')
#ylab = 'Probability of approach'
#dat <- metrics$approach
#make_plot(dat, plotpath, ylab, ylim = c(0.4,0.6),abline_y = 0.5)

#TODO:
#look into autocorrelation, consider downsampling 
#model with autocorrelation model to get error bars 
