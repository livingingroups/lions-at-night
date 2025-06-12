#What is lion cohesion + coordination like during the night?
#What are the relevant spatial scales associated with lion coordination?
#This script so far has 3 parts:
#Part 1 looks at spatial scales of coordination by plotting various dyadic 
#coordination-related metrics (e.g. heading correlation, change in dyadic distance, etc.)
#as a function of distance between pairs.
#Part 2 makes some more basic descriptive plots of speed, spread, polarization, etc. 
#and their relationship with time of day.
#Part 3 looks into leader/follower relationships during the night. This is defined based
#on the probability that one individual (the follower) moves in the direction of another's
#current location (the leader)

library(cocomo)
library(lubridate)
library(fields)

#PARAMETERS - MODIFY THESE

#directories
dir <- '~/Dropbox/lions_at_night/' #path to project directory
speed_dt <- 10 #time difference to use for computing speeds and headings (in units of timesteps)
min_speed_to_compute_heading <- .5 #minimum speed needed to have a heading (otherwise headings are set to NA)
start_times <- as.POSIXct(paste(seq.Date(date('2023-05-06'),date('2023-06-15'),by=1), '16:00:00'), tz = 'UTC') #timestamps of times to start each analysis period for spatial scales analyses
end_times <- as.POSIXct(paste(seq.Date(date('2023-05-07'),date('2023-06-16'),by=1), '06:00:00'), tz = 'UTC') #timestamps of times to end each analysis period for spatial scales analyses
dist_bins <- c(0, 10^seq(0.6,4,.4)) #distance between dyads bins to use for spatial scales analyses
fit_gams <- F
min_points_to_plot <- 3*60*60 #minimum of 3 hrs of data needed to plot things

#------------------------------------------------
#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#FUNCTIONS

#helper function to get relevant statistics of the distribution of values for a set of bins (binned by another variable x)
#vals = values for which you are getting the distribution of (e.g. heading correlation)
#x = what you are binning on (e.g. dyadic distance)
#bins = bins to use for binning x
#the dimensions of vals and x must be the same
get_distrib_statistics <- function(vals, x, bins, fit_gam = F, gam_family = 'gaussian', gam_k = 5){
  out <- list()
  out$bins <- bins
  out$mean <- out$median <- out$q75 <- out$q25 <- out$q025 <- out$q975 <- rep(NA, length(bins)-1)
  #loop over all bins
  for(i in 1:(length(bins)-1)){
    
    #get indexes to data in that bin
    idxs <- which(x >= bins[i] & x < bins[i+1])
    
    #get relevant stats
    out$mean[i] <- mean(vals[idxs], na.rm=T)
    out$median[i] <- median(vals[idxs], na.rm=T)
    out$q75[i] <- quantile(vals[idxs], 0.75, na.rm=T)
    out$q25[i] <- quantile(vals[idxs], 0.25, na.rm=T)
    out$q025[i] <- quantile(vals[idxs], 0.025, na.rm=T)
    out$q975[i] <- quantile(vals[idxs], 0.975, na.rm=T)
  }
  
  if(fit_gam){
    mod <- gam(formula = val ~ s(x, 3), family = gam_family, data = data.frame(val = c(vals), x = c(x)))
    out$gam <- mod
  }
  
  return(out)
}

#helper function to make plots of different metrics
make_plot <- function(dat, plotpath, ylab, ylim = NULL, xlab = 'Distance apart (m)', abline_y = NULL, logx = T, plot_means = T, plot_medians = T, plot_IQR = T){
  if(is.null(ylim)){
    ylim <- c(min(dat$q25,na.rm=T), max(dat$q75, na.rm=T))
  }
  
  png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
  mids <- (dat$bins[2:length(dat$bins)] + dat$bins[1:(length(dat$bins)-1)]) / 2
  if(logx){
    plot(mids, dat$mean, xlab = xlab, ylab = ylab, pch = 19, col = '#00000000', cex = 1.5, ylim = ylim ,log='x')
  } else{
    plot(mids, dat$mean, xlab = xlab, ylab = ylab, pch = 19, col = '#00000000', cex = 1.5, ylim = ylim)
    
  }
  if(!is.null(abline_y)){
    abline(h=abline_y, lty = 2)
  }
  
  if(!is.null(dat$gam)){
    newdat <- data.frame(x = seq(min(dat$bins),max(dat$bins), length.out=10000))
    preds <- predict(dat$gam, newdata = newdat, type = 'response')
    lines(newdat$x, preds, lwd = 3)
  }
  
  if(plot_IQR){
    arrows(mids, dat$q25, mids, dat$q75, lwd = 2, code = 3, length = 0.1, angle = 90)
  }
  if(plot_medians){
    points(mids, dat$median, pch = 19, col = 'black',cex = 1.5)
  }
  if(plot_means){
    points(mids, dat$mean, pch = 19, col = 'gray', cex = 1.5)
  }
  dev.off()
}

#----LOAD AND PROCESS DATA----

#Load data
gpsfile <- paste0(dir, 'data_namibia/processed/lion_omukutu_xy_level1.RData')
load(gpsfile) 

#basic info
n_inds <- nrow(xs)
n_times <- ncol(xs)

#Get periods of data to use for spatial scales analyses
periods <- data.frame(start_time = start_times, end_time = end_times)
periods$t0_idx <- match(periods$start_time, timestamps)
periods$tf_idx <- match(periods$end_time, timestamps)

#relevant time indexes for spatial scales analyses
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
  
  #remove heading when speed is too low
  heads[i,which(speeds[i,] < min_speed_to_compute_heading)] <- NA
}

#dyadic distances
dyad_dists <- cocomo::get_group_dyadic_distances(xs, ys)

#group polarization (computed for all instances where at least 5 individuals were tracked)
polarization <- cocomo::get_group_polarization(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)

#group speed (computed for all instances where at least 5 individuals were tracked)
out <- cocomo::get_group_heading_and_speed(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)
group_speed <- out$speeds

#get mean dyadic distance (computed for all instances where at least 5 individuals were tracked)
mean_dyad_dist <- apply(dyad_dists, 3, FUN = mean, na.rm=T)
n_tracked <- colSums(!is.na(xs))
mean_dyad_dist[which(n_tracked < n_inds - 1)] <- NA

#dyadic distance changes, speed differences, and heading correlations between pairs of individuals
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

#heading toward (or away from) conspecifics' previous locations
head_twd <- array(NA, dim = c(n_inds, n_inds, n_times))
for(i in 1:n_inds){
  for(j in 1:n_inds){
    if(i != j){
      
      #heading of i
      head_i <- heads[i,] 
      
      #vector pointing from i to j
      dx_ij <- (xs[j,] - xs[i,]) / dyad_dists[i,j,]
      dy_ij <- (ys[j,] - ys[i,]) / dyad_dists[i,j,]
      
      #angle between these vectors
      head_twd[i,j,] <- acos(cos(head_i)*dx_ij + sin(head_i)*dy_ij)
    }
  }
}


#get night data only (set non-night indexes to NA)
dyad_dists_night <- dyad_dists
head_corrs_night <- head_corrs
speed_diffs_night <- speed_diffs
log_speed_diffs_night <- log_speed_diffs
dyad_dist_changes_night <- dyad_dist_changes
head_twd_night <- head_twd

dyad_dists_night[,,idxs_to_remove] <- NA
head_corrs_night[,,idxs_to_remove] <- NA
speed_diffs_night[,,idxs_to_remove] <- NA
log_speed_diffs_night[,,idxs_to_remove] <- NA
dyad_dist_changes_night[,,idxs_to_remove] <- NA
head_twd_night[,,idxs_to_remove] <- NA

#remove data when neither individual is moving from dyadic distance changes
dyad_dist_changes_when_moving_night <- dyad_dist_changes_night
for(i in 1:(n_inds-1)){
  for(j in 2:n_inds){
    slow_i <- speeds[i,] < min_speed_to_compute_heading
    slow_j <- speeds[j,] < min_speed_to_compute_heading
    idxs_slow_both <- which(slow_i & slow_j)
    dyad_dist_changes_when_moving_night[i,j,idxs_slow_both] <- NA
    dyad_dist_changes_when_moving_night[j,i,idxs_slow_both] <- NA
  }
}

#calculate mean, median, iqr of metrics during relevant periods only, as a function of distance apart
metrics <- list()
metrics$dist_bins <- dist_bins

#get statistics of the distribution for each metric within each bin
metrics$head_corr <- get_distrib_statistics(head_corrs_night, dyad_dists, dist_bins)
metrics$speed_diff <- get_distrib_statistics(speed_diffs_night, dyad_dists, dist_bins)
metrics$log_speed_diff <- get_distrib_statistics(log_speed_diffs_night, dyad_dists, dist_bins)
metrics$dyad_dist_change <- get_distrib_statistics(dyad_dist_changes_night, dyad_dists, dist_bins)
metrics$dyad_dist_change_when_moving <- get_distrib_statistics(dyad_dist_changes_when_moving_night, dyad_dists, dist_bins)
metrics$approach <- get_distrib_statistics(dyad_dist_changes_night[which(dyad_dist_changes!=0)] < 0, dyad_dists, dist_bins)
metrics$approach_when_moving <- get_distrib_statistics(dyad_dist_changes_when_moving_night < 0, dyad_dists, dist_bins)
metrics$ang_between_heads <- get_distrib_statistics(acos(head_corrs_night)*180/pi, dyad_dists, dist_bins)
metrics$head_twd <- get_distrib_statistics(head_twd_night, dyad_dists, dist_bins)

#-----PLOTS-----

#----Part 1: Spatial scales analyses----

#Plot: Angle between headings 
plotpath <- paste0(dir, 'plots/spatial_scales/ang_between_heads.png')
ylab <- "Angle between headings (degrees)"
ylim <- c(0,180)
dat <- metrics$ang_between_heads
make_plot(dat, plotpath, ylab, ylim, abline_y = 90)

#Plot: Heading correlation
plotpath <- paste0(dir, 'plots/spatial_scales/heading_corr.png')
ylab <- "Heading correlation"
ylim <- c(-1,1)
dat <- metrics$head_corr
make_plot(dat, plotpath, ylab, ylim, abline_y = 0)

#Plot: Mean difference in speed as a function of distance 
plotpath <- paste0(dir, 'plots/spatial_scales/abs_speed_diff.png')
ylab = 'Absolute speed difference (m/s)'
dat <- metrics$speed_diff
make_plot(dat, plotpath, ylab)

#Plot: Change in dyadic distance
plotpath <- paste0(dir, 'plots/spatial_scales/dyad_dist_change.png')
ylab = 'Change in dyadic distance (m)'
dat <- metrics$dyad_dist_change
make_plot(dat, plotpath, ylab, abline_y = 0)

#Plot: Change in dyadic distance, only when moving
plotpath <- paste0(dir, 'plots/spatial_scales/dyad_dist_change_when_moving.png')
ylab = 'Change in dyadic distance (m)'
dat <- metrics$dyad_dist_change_when_moving
make_plot(dat, plotpath, ylab, abline_y = 0)


#Plot: Probability of approaching (when moving)
plotpath <- paste0(dir, 'plots/spatial_scales/p_approach_when_moving.png')
ylab = 'Probability of approach'
dat <- metrics$approach_when_moving
make_plot(dat, plotpath, ylab, ylim = c(0,0.6),abline_y = 0.5, plot_medians = F, plot_IQR = F)

#Plot: Heading towards
plotpath <- paste0(dir, 'plots/spatial_scales/head_towards_conspecific.png')
ylab = 'Heading relative to vector toward conspecific'
dat <- metrics$head_twd
make_plot(dat, plotpath, ylab, ylim = NULL,abline_y = pi/2, plot_medians = T, plot_IQR = T)

#----Part 2: Basic info during the night----

hour_of_day <- lubridate::hour(timestamps)

#Get metrics by hour
metrics_by_hour <- list()
metrics_by_hour$mean_dyad_dist <- get_distrib_statistics(mean_dyad_dist, hour_of_day, seq(-.5,24,1))
metrics_by_hour$polarization <- get_distrib_statistics(polarization, hour_of_day, seq(-.5,24,1))

#Plot: Distribution of log(dyadic distances) between individuals
png(filename = paste0(dir, 'plots/overall_hists/hist_log_dyad_dist.png'), width = 8, height = 6, units = 'in', res = 300)
histo <- hist(log(dyad_dists, 10), plot=T, breaks=80, xlab = 'Dyadic distance (m) - log bins',main='', xlim = c(0,4), freq = F, xaxt ='n')
axis(1, at = seq(0,4,1), labels = 10^seq(0,4,1))
dev.off()

#Plot: Mean dyadic distance by hour
plotpath <- paste0(dir, 'plots/hourly_metrics/mean_dyad_dist_by_hr.png')
dat <- metrics_by_hour$mean_dyad_dist
make_plot(dat, plotpath, ylab = 'Mean dyadic distance (m)', logx=F, xlab = 'Hour UTC', plot_means = F)

#Plot: Polarization by hour
plotpath <- paste0(dir, 'plots/hourly_metrics/pol_by_hr.png')
dat <- metrics_by_hour$polarization
make_plot(dat, plotpath, ylab = 'Group polarization', logx=F, xlab = 'Hour UTC', plot_means = F, ylim =c(0,1))

#Plot: histogram of group polarization
png(filename = paste0(dir, 'plots/overall_hists/hist_polarization.png'), width = 8, height = 6, units = 'in', res = 300)
histo <- hist(polarization, plot=T, breaks=seq(0,1,.01), xlab = 'Group polarization',main='', xlim = c(0,1), freq = F)
dev.off()

#Plot: cumulative histogram of group speed and individual speed
png(filename = paste0(dir, 'plots/overall_hists/hist_speed.png'), width = 8, height = 6, units = 'in', res = 300)
speedx <- seq(0,max(group_speed,na.rm=T)+1,.01)
histo <- hist(group_speed, breaks = speedx, plot = F)
cumhist <- cumsum(histo$counts) / sum(histo$counts)
plot(histo$mids, cumhist,log='x', xlim=c(.1,10), type = 'l', lwd = 3, xlab = 'Speed (m/s)', ylab = 'Cumulative probability', ylim=c(0.01,1))
speedx <- seq(0,max(speeds,na.rm=T)+1,.01)
histo <- hist(speeds, breaks = speedx, plot = F)
cumhist <- cumsum(histo$counts) / sum(histo$counts)
lines(histo$mids, cumhist, col = 'red', lwd = 3)
abline(v=1.3, lty = 2)
text(1.1, 0.5, 'Person walking', srt = 90)
abline(v=5, lty = 2)
text(4.3, .5, 'Person running', srt = 90)
legend('bottomleft', legend = c('Individual','Group centroid'), col = c('red','black'),lwd=c(3,3))
dev.off()

#Plot: Speed distribution as a function of time of day (individual)
png(filename = paste0(dir, 'plots/hourly_metrics/speed_distrib_by_hr.png'), width = 8, height = 6, units = 'in', res = 300)
speedx <- seq(0,max(speeds,na.rm=T)+1,.001)
plot(NULL, xlim = c(.1,10), ylim = c(0,1), log = 'x', xlab = 'Individual speed (m/s)', ylab = 'Cumulative probability')
hour_bins <- seq(0,24,2)
cols <- viridis(length(hour_bins)-1)
for(hour in 1:(length(hour_bins)-1)){
  idxs <- which(hour_of_day >= hour_bins[hour] & hour_of_day < hour_bins[hour+1])
  histo <- hist(speeds[,idxs], breaks = speedx, plot = F)
  cumhist <- cumsum(histo$counts) / sum(histo$counts)
  lines(histo$mids, cumhist, col = cols[hour], lwd = 2)
}
labs <- paste0(hour_bins[1:(length(hour_bins)-1)],'-', hour_bins[2:length(hour_bins)],' UTC')
legend('bottomright', legend = labs, col = cols, lwd = 2)
dev.off()

#Plot: Speed distribution as a function of time of day (group)
png(filename = paste0(dir, 'plots/hourly_metrics/centr_speed_distrib_by_hr.png'), width = 8, height = 6, units = 'in', res = 300)
speedx <- seq(0,max(speeds,na.rm=T)+1,.001)
plot(NULL, xlim = c(.1,10), ylim = c(0,1), log = 'x', xlab = 'Centroid speed (m/s)', ylab = 'Cumulative probability')
hour_bins <- seq(0,24,2)
cols <- viridis(length(hour_bins)-1)
for(hour in 1:(length(hour_bins)-1)){
  idxs <- which(hour_of_day >= hour_bins[hour] & hour_of_day < hour_bins[hour+1])
  histo <- hist(group_speed[idxs], breaks = speedx, plot = F)
  cumhist <- cumsum(histo$counts) / sum(histo$counts)
  lines(histo$mids, cumhist, col = cols[hour], lwd = 2)
}
labs <- paste0(hour_bins[1:(length(hour_bins)-1)],'-', hour_bins[2:length(hour_bins)],' UTC')
legend('bottomright', legend = labs, col = cols, lwd = 2)
dev.off()

#---Part 3: Leader/follower analyses---

#Plot: Pairwise leader/follower relations
#for each pair, get probability of i moving toward j
p_toward_ij <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    if(i!=j){
      p_toward_ij[i,j] <- mean(head_twd_night[i,j,] < pi/2, na.rm=T)
    }
  }
}

plotpath <- paste0(dir, 'plots/dyadic/prob_heading_towards.png')
png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
image.plot(p_toward_ij, col = viridis(256), xlab = 'Follower', ylab = 'Leader', xaxt ='n',yaxt='n')
axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code)
axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code)
dev.off()

#Plots: Heading toward as a function of initial distance apart
#Approach as a function of distance (x axis), multi-panel by leader
p_approach_vs_dist <- uppers <- lowers <- tots <- array(NA, dim = c(n_inds, n_inds, length(dist_bins)-1))
for(d in 1:(length(dist_bins)-1)){
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      if(i != j){
        curr_twd <- head_twd_night[i,j,] < pi/2
        curr_dists <- dyad_dists[i,j,]
        idxs <- which(curr_dists >= dist_bins[d] & curr_dists < dist_bins[d+1])
        p_approach_vs_dist[i,j,d] <- mean(curr_twd[idxs], na.rm=T)
        n_approach <- sum(curr_twd[idxs]==1, na.rm=T)
        n_tot <- sum(!is.na(curr_twd[idxs]))
        if(n_tot > 0){
          lowers[i,j,d] <- binom.test(n_approach, n_tot)$conf.int[1]
          uppers[i,j,d] <- binom.test(n_approach, n_tot)$conf.int[2]
        }
        tots[i,j,d] <- n_tot
      }
    }
  }
}

#by leader
plotpath <- paste0(dir, 'plots/dyadic/prob_approach_vs_dist_by_leader.png')
png(filename = plotpath, width = 12, height = 6, units = 'in', res = 300)
par(mfrow = c(2,3), mar = c(3,3,1,1))
mids <- dist_bins[1:(length(dist_bins)-1)] + diff(dist_bins)/2
cols <- viridis(n_inds)
cols_err <- paste0(substr(cols, 1,7),'33')
for(j in 1:n_inds){
  plot(NULL, xlim = c(1, dist_bins[length(dist_bins)-1]), ylim = c(0,1), log = 'x', main = paste0('Leader = ',ids$code[j]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
  for(i in 1:n_inds){
    if(i!=j){
      curr <- p_approach_vs_dist[i,j,]
      curr[which(tots[i,j,]<min_points_to_plot)] <- NA
      lines(mids, curr, lwd = 2, col = cols[i])
    }
  }
  abline(h=0.5, lty = 2)
  legend('bottomleft', col = cols, lwd = 2, legend = paste0('Follower = ',ids$code))
}
dev.off()

#by follower
plotpath <- paste0(dir, 'plots/dyadic/prob_approach_vs_dist_by_follower.png')
png(filename = plotpath, width = 12, height = 6, units = 'in', res = 300)
par(mfrow = c(2,3), mar = c(3,3,1,1))
mids <- dist_bins[1:(length(dist_bins)-1)] + diff(dist_bins)/2
cols <- viridis(n_inds)
cols_err <- paste0(substr(cols, 1,7),'33')
for(i in 1:n_inds){
  plot(NULL, xlim = c(1, dist_bins[length(dist_bins)-1]), ylim = c(0,1), log = 'x', main = paste0('Follower = ',ids$code[i]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
  for(j in 1:n_inds){
    if(i!=j){
      curr <- p_approach_vs_dist[i,j,]
      curr[which(tots[i,j,]<min_points_to_plot)] <- NA
      lines(mids, curr, lwd = 2, col = cols[j])
    }
  }
  abline(h=0.5, lty = 2)
  legend('bottomleft', col = cols, lwd = 2, legend = paste0('Leader = ',ids$code))
}
dev.off()

#Plot: Approach probability for all dyads (heatmap), multi-panel by distance bin
plotpath <- paste0(dir, 'plots/dyadic/approach_heatmaps_by_dist.png')
png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
par(mfrow = c(3,3), mar = c(3,3,1,1))
for(d in 1:(length(dist_bins)-1)){
  curr <- p_approach_vs_dist[,,d]
  curr[which(tots[,,d]<min_points_to_plot)] <- NA
  image.plot(curr, xaxt='n', yaxt='n' ,zlim=c(0,1),col=viridis(256),xlab='Follower',ylab='Leader', main = paste0(round(dist_bins[d]),'-',round(dist_bins[d+1]),' m'))
  axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code, las = 2,cex.axis=0.5)
  axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code, las=2,cex.axis=0.5)
  
}
dev.off()

#Plot: Approach probability for the most relevant "middle" distance window - 10 m to 500 m.
p_approach_middle <- n_middle <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    idxs <- which(dyad_dists[i,j,] >= 10 & dyad_dists[i,j,] < 500)
    p_approach_middle[i,j] <- mean(head_twd_night[i,j,idxs] < pi/2, na.rm=T)
    n_middle[i,j] <- sum(!is.na(head_twd_night[i,j,idxs]))
  }
}

plotpath <- paste0(dir, 'plots/dyadic/approach_heatmap_10-500m.png')
png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
image.plot(p_approach_middle, xaxt='n', yaxt='n' ,zlim=c(0,1),col=viridis(256),xlab='Follower',ylab='Leader', main = 'P(approach) when 10-500 m apart')
axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code, las = 2)
axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code, las=2)
y <- c(matrix(rep(seq(0,1,length.out=n_inds), each = n_inds), nrow = n_inds, ncol = n_inds))
x <- c(matrix(rep(seq(0,1,length.out=n_inds), n_inds), nrow = n_inds, ncol = n_inds))
text(x,y,paste0(round(c(p_approach_middle)*100),'%'))
dev.off()

#Plot: correlation between leading vs following rates
plotpath <- paste0(dir, 'plots/dyadic/leader_follower_correlation_10-500m.png')
png(filename = plotpath, width = 6, height = 6, units = 'in', res = 300)
plot(c(p_approach_middle),c(t(p_approach_middle)), cex = 1, pch = 19, xlab = 'Probability of approaching',ylab = 'Probability of being approached', xlim = c(0,1),ylim=c(0,1), asp = 1)
lines(c(0,1),c(1,0),lty =2)
dev.off()

#TODO:
#look into autocorrelation, consider downsampling 
#model with autocorrelation model to get error bars 
#nearest neighbor
#different time scales - 5 sec, 10 sec, 30, sec, 1 min ,5, 10
#remove dyad dist changes if both below a speed threshold

data <- data.frame(dist = c(dyad_dists), approach = c(dyad_dist_changes_when_moving_night < 0))
data <- data[complete.cases(data),]
data <- data[which(data$dist < 10000),]
test <- gam(approach ~ s(dist,3), data = data, family = 'binomial')
newx <- seq(1,5000)
predicted <- predict(test, newdata = data.frame(dist = newx), type = 'response',se=T)
plot(newx, predicted, type = 'l', lwd = 3, log = 'x')
