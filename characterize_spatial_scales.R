#What is lion cohesion + coordination like during the night?
#What are the relevant spatial scales associated with lion coordination?
#How do individuals vary in their cohesive behavior?
#This script so far runs 4 sets of analyses:
#Analysis 1 looks at spatial scales of coordination by plotting various dyadic 
#coordination-related metrics (e.g. heading correlation, change in dyadic distance, etc.)
#as a function of distance between pairs, across all dyads.
#Analysis 2 makes some more basic descriptive plots of speed, spread, polarization, etc. 
#and their relationship with time of day.
#Analysis 3 looks at cohesion on a dyadic level (i.e. thinking about relationships between specific dyads,
#so a kind of leader-follower type analysis).
#This is defined based on the probability that one individual (the "ego") moves in the direction of another's
#(the "other") current location. In some places in the script I call the "ego" the "follower" and the
#"other" the "leader". I couldn't be bothered to fix it, sorry.
#Analysis 4 is similar to Analysis 3, but we now also consider what the "other" is doing, i.e. whether 
#it is currently stationary, moving toward the ego, or moving away from the ego.
#There is more description of all these analyses in the lion project notebook

#The script is organized such that the first part (LOAD AND PROCESS DATA) does the computations and the second part (PLOTS) 
#is basically just plotting (though there are a few computations in there because I didn't manage to clean this up fully)
#From the first part, the data frame approach_data is optionally saved (this is for the analyses in Analysis 4)

#-------------LIBRARIES--------------
library(cocomo)
library(lubridate)
library(fields)
library(ggplot2)
library(viridis)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#------------PARAMETERS-----------

dir <- '~/Dropbox/lions_at_night/' #path to project directory
speed_dt <- 10 #time difference to use for computing speeds and headings (in units of timesteps i.e. seconds)
min_speed_to_compute_heading <- .5 #minimum speed needed to have a heading (otherwise headings are set to NA). Also, when individuals are below this speed they are considered 'stationary' for analysis 4.
dist_bins <- c(0, 10^seq(0.6,4,.4)) #distance between dyads bins to use for spatial scales analyses
min_points_to_plot <- 3*60*60 #minimum of 3 hrs of data needed to plot things (currently only used in analysis 3 and 4)
fit_gams <- F #this didn't work well, just leave it as FALSE

start_times <- as.POSIXct(paste(seq.Date(date('2023-05-06'),date('2023-06-15'),by=1), '16:00:00'), tz = 'UTC') #timestamps of times to start each analysis period (night) - replace eventually with info on when the group left its rest spot
end_times <- as.POSIXct(paste(seq.Date(date('2023-05-07'),date('2023-06-16'),by=1), '06:00:00'), tz = 'UTC') #timestamps of times to end each analysis period (night) - replace eventually with info on when the group left its rest spot

savename_df <- 'data_namibia/processed/cohesion_data.RData' #where to save the data frame 'approach_data' (from analysis 4) - if NULL, nothing will be saved
make_plots <-  T #whether to generate plots or not

#-------------FUNCTIONS----------------

#helper function to get relevant statistics of the distribution of values for a set of bins (binned by another variable x)
#vals = values for which you are getting the distribution of (e.g. heading correlation)
#x = what you are binning on (e.g. dyadic distance)
#bins = bins to use for binning x
#the dimensions of vals and x must be the same
get_distrib_statistics <- function(vals, x, bins, fit_gam = F, gam_family = 'gaussian', gam_k = 5){
  out <- list()
  out$bins <- bins
  out$mean <- out$median <- out$q75 <- out$q25 <- out$q025 <- out$q975 <- out$n <- rep(NA, length(bins)-1)
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
    out$n[i] <- sum(!is.na(vals[idxs]))
  }
  
  if(fit_gam){
    mod <- gam(formula = val ~ s(x, 3), family = gam_family, data = data.frame(val = c(vals), x = c(x)))
    out$gam <- mod
  }
  
  return(out)
}

#helper function to make plots of different metrics vs distance
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

#Plotting function for 'interactions' (i.e. following, converging, joining)
#other_behav should be either NULL (use all data), 'toward', 'away', or 'stationary'
#stratify_by gives whether to stratify by dyad, can be 'ego','other' or NULL (if you want to make the plot aggregating across all dyads)
make_approach_vs_dist_plot <- function(approach_data, dist_bins, plotpath, other_behav = NULL, include_moving_only = F, stratify_by = NULL, min_points_to_plot = 60*60){
  
  #if indexes to subset by aren't specified, use all rows of data
  if(is.null(other_behav)){
    idxs_to_use <- 1:nrow(approach_data)
  } else{
    idxs_to_use <- which(approach_data$past_behav_j == other_behav)
  }
  
  #subset to only data we need to use
  dat <- approach_data[idxs_to_use,]
  
  #get midpoints of bins
  mids <- dist_bins[1:(length(dist_bins)-1)] + diff(dist_bins)/2
  
  #main label
  if(is.null(other_behav)){
    main_lab <- 'All data'
  } else{
    if(other_behav == 'toward'){
      main_lab <- 'Converging'
    }
    if(other_behav == 'away'){
      main_lab = 'Following'
    }
    if(other_behav == 'stationary'){
      main_lab = 'Joining'
    }
  }
  
  if(is.null(stratify_by)){
    #get means by distance bin
    if(include_moving_only){
      toward <- get_distrib_statistics(vals = dat$i_approaches_j_given_moving, x = dat$dyad_dist, bins = dist_bins)
    } else{
      toward <- get_distrib_statistics(vals = dat$i_approaches_j, x = dat$dyad_dist, bins = dist_bins)
    }
    
    #make plot
    png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
    if(include_moving_only){
      plot(mids, toward$mean, pch = 19, cex = 2, col = 'gray', log = 'x', xlab = 'Distance apart (m)', ylab = 'Probability of approach (given moved)', main = main_lab)
      abline(h=0.5, lty = 2)
    } else{
      plot(mids, toward$mean, pch = 19, cex = 2, col = 'gray', log = 'x', xlab = 'Distance apart (m)', ylab = 'Probability of approach', main = main_lab)
    }
    dev.off()
  } else{
    n_inds <- max(approach_data$i_ego)
    n_bins <- length(dist_bins)-1
    approach_probs <- ns <- array(NA, dim = c(n_inds, n_inds, n_bins))
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        if(i!=j){
          idxs <- which(dat$i_ego==i & dat$j_other==j)
          #get means by distance bin
          if(include_moving_only){
            toward <- get_distrib_statistics(vals = dat$i_approaches_j_given_moving[idxs], x = dat$dyad_dist[idxs], bins = dist_bins)
          } else{
            toward <- get_distrib_statistics(vals = dat$i_approaches_j[idxs], x = dat$dyad_dist[idxs], bins = dist_bins)
          }
          approach_probs[i,j,] <- toward$mean
          ns[i,j,] <- toward$n
        }
      }
    }
    approach_probs[which(ns < min_points_to_plot)] <- NA
    ylims <- c(min(approach_probs,na.rm=T), max(approach_probs,na.rm=T))
    
    #make the plot
    png(filename = plotpath, width = 12, height = 6, units = 'in', res = 300)
    par(mfrow = c(2,3), mar = c(3,3,1,1))
    cols <- viridis(n_inds)
    
    if(stratify_by == 'other'){
      for(j in 1:n_inds){
        plot(NULL, xlim = c(1, dist_bins[length(dist_bins)-1]), ylim = ylims, log = 'x', main = paste0('Other = ',ids$code[j]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
        for(i in 1:n_inds){
          if(i!=j){
            curr <- approach_probs[i,j,]
            lines(mids, curr, lwd = 2, col = cols[i])
          }
        }
        if(include_moving_only){
          abline(h=0.5, lty = 2)
        }
        legend('bottomleft', col = cols, lwd = 2, legend = paste0('Ego = ',ids$code))
      }
    }
    if(stratify_by == 'ego'){
      for(i in 1:n_inds){
        plot(NULL, xlim = c(1, dist_bins[length(dist_bins)-1]), ylim = ylims, log = 'x', main = paste0('Ego = ',ids$code[i]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
        for(j in 1:n_inds){
          if(i!=j){
            curr <- approach_probs[i,j,]
            lines(mids, curr, lwd = 2, col = cols[j])
          }
        }
        if(include_moving_only){
          abline(h=0.5, lty = 2)
        }
        legend('bottomleft', col = cols, lwd = 2, legend = paste0('Other = ',ids$code))
      }
    }
    dev.off()
    
    
    
  }
  
}


#----LOAD AND PROCESS DATA----

print('loading and processing data')

#Load data
gpsfile <- paste0(dir, 'data_namibia/processed/lion_omukutu_xy_level1.RData')
load(gpsfile) 

#basic info
n_inds <- nrow(xs)
n_times <- ncol(xs)

#Get periods of data to use for spatial scales analyses
#To Do: replace with stop/start periods defined by Gen's analysis
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

#headings and speeds (future)
heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
for(i in 1:n_inds){
  out <- cocomo::get_heading_and_speed_temporal(x_i = xs[i, ], y_i = ys[i,], t_window = speed_dt, forward = T)
  heads[i,] <- out$heads
  speeds[i,] <- out$speeds
  
  #remove heading when speed is too low
  heads[i,which(speeds[i,] < min_speed_to_compute_heading)] <- NA
}

#past headings and speeds
heads_past <- matrix(NA, nrow = n_inds, ncol = n_times)
speeds_past <- matrix(NA, nrow = n_inds, ncol = n_times)
for(i in 1:n_inds){
  out <- cocomo::get_heading_and_speed_temporal(x_i = xs[i,], y_i = ys[i,], t_window = speed_dt, forward = F)
  heads_past[i,] <- out$heads
  speeds_past[i,] <- out$speeds
  
  #remove heading when speed is too low
  heads_past[i,which(speeds_past[i,] < min_speed_to_compute_heading)] <- NA
}

#dyadic distances
dyad_dists <- cocomo::get_group_dyadic_distances(xs, ys)

#group polarization (computed for all instances where at least 5 individuals were tracked, otherwise NA)
polarization <- cocomo::get_group_polarization(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)

#group speed (computed for all instances where at least 5 individuals were tracked, otherwise NA)
out <- cocomo::get_group_heading_and_speed(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)
group_speed <- out$speeds

#get mean dyadic distance (computed for all instances where at least 5 individuals were tracked, otherwise NA)
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

#heading relative to others' previous locations
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

#get night data only for some plots (set non-night indexes to NA)
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

#----Data frame to hold all 'interactions' data at night----
#i is the "ego", i.e. the one whose behavior in response to "other" we are considering
#j is the "other, i.e. the reference individual relative to whom i's behavior is being considered
approach_data <- data.frame(i_ego = rep(rep(1:n_inds, n_inds),n_times), 
                            j_other = rep(rep(1:n_inds, each = n_inds),n_times),
                            tidx = rep(1:n_times, each = n_inds*n_inds))

#remove non-night columns (otherwise we might run out of memory :/)
approach_data <- approach_data[which(approach_data$tidx %in% relevant_t_idxs),]

#columns in the approach_data data frame:
#x_i = x location of individual i
#y_i = y location of individual i
#x_j = x location of individual j
#y_j = y location of individual j
#head_i_fut - future heading of i (direction vector pointing from my position now to my position speed_dt seconds into the future)
#head_i_past - past heading of i (direction vector pointing from my position in the past (speed_dt seconds ago) to my position now)
#head_j_fut <- future heading of j
#head_j_past <- past heading of j
#speed_i_fut <- speed of i in the future
#speed_j_fut <- speed of j in the future
#speed_i_past <- speed of i in the past
#speed_j_past <- speed of j in the past
#dyad_dist - distance between i and j
approach_data$x_i <- xs[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$y_i <- ys[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$x_j <- xs[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$y_j <- ys[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$head_i_fut <- heads[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$head_j_fut <- heads[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$head_i_past <- heads_past[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$head_j_past <- heads_past[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$speed_i_fut <- speeds[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$speed_j_fut <- speeds[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$speed_i_past <- speeds_past[cbind(approach_data$i_ego, approach_data$tidx)]
approach_data$speed_j_past <- speeds_past[cbind(approach_data$j_other, approach_data$tidx)]
approach_data$dyad_dist <- dyad_dists[cbind(approach_data$i_ego, approach_data$j_other, approach_data$tidx)]

#columns in the approach_data data frame:
#unit vector pointing from i to j and from j to i (for use below)
dx_itoj <- (approach_data$x_j - approach_data$x_i) / approach_data$dyad_dist
dy_itoj <- (approach_data$y_j - approach_data$y_i) / approach_data$dyad_dist
dx_jtoi <- (approach_data$x_i - approach_data$x_j) / approach_data$dyad_dist
dy_jtoi <- (approach_data$y_i - approach_data$y_j) / approach_data$dyad_dist

#columns in the approach_data data frame:
#angle between past heading of j_other and vector pointing from j_other to i_ego (0 means j was moving directly toward i's position, pi means j was moving directly away)
approach_data$j_past_head_twd_i <- acos(cos(approach_data$head_j_past)*dx_jtoi + sin(approach_data$head_j_past)*dy_jtoi)

#columns in the approach_data data frame:
#angle between future heading of i_ego and vector pointing from i_ego to j_other (0 means i_ego moves directly toward j's position, pi means i moves directly away)
approach_data$i_fut_head_twd_j <- acos(cos(approach_data$head_i_fut)*dx_itoj + sin(approach_data$head_i_fut)*dy_itoj)

#columns in the approach_data data frame:
#past behavior of j_other in categories - 'away' from i_ego, 'toward' i_ego or 'stationary' (not moving more than min_speed_to_compute_heading)
approach_data$past_behav_j <- NA
approach_data$past_behav_j[which(approach_data$j_past_head_twd_i >= pi/2 & approach_data$speed_j_past >= min_speed_to_compute_heading)] <- 'away'
approach_data$past_behav_j[which(approach_data$j_past_head_twd_i < pi/2 & approach_data$speed_j_past >= min_speed_to_compute_heading)] <- 'toward'
approach_data$past_behav_j[which(approach_data$speed_j_past < min_speed_to_compute_heading)] <- 'stationary'

#columns in the approach_data data frame:
#future behavior of i_ego - 'away' from j_other, 'toward' j_other, or 'statinoary' (not moving more than min_speed_to_compute_heading)
approach_data$fut_behav_i <- NA
approach_data$fut_behav_i[which(approach_data$i_fut_head_twd_j >= pi/2 & approach_data$speed_i_fut >= min_speed_to_compute_heading)] <- 'away'
approach_data$fut_behav_i[which(approach_data$i_fut_head_twd_j < pi/2 & approach_data$speed_i_fut >= min_speed_to_compute_heading)] <- 'toward'
approach_data$fut_behav_i[which(approach_data$speed_i_fut < min_speed_to_compute_heading)] <- 'stationary'

#columns in the approach_data data frame:
#future behavior of i_ego as a binary variable, including stationary as 'not' doing the behavior
#this will be 1 if i_ego approaches j_other and 0 if they either remain stationary or move away
approach_data$i_approaches_j <- NA
approach_data$i_approaches_j[which(!is.na(approach_data$fut_behav_i))] <- 0
approach_data$i_approaches_j[which(approach_data$fut_behav_i == 'toward')] <- 1

#columns in the approach_data data frame:
#future behavior of i_ego as a binary variable, excluding when they are stationary (i.e. stationary is considered NA)
#this will be 1 if i_ego approaches j_other, 0 if they move away, and NA if they are stationary
approach_data$i_approaches_j_given_moving <- NA
approach_data$i_approaches_j_given_moving[which(approach_data$fut_behav_i == 'away')] <- 0
approach_data$i_approaches_j_given_moving[which(approach_data$fut_behav_i == 'toward')] <- 1

#save data frame
if(!is.null(savename_df)){
  print('saving cohesion dataframe (approach_data)')
  save(list = c('approach_data'), file = paste0(dir,savename_df))
}

#-----PLOTS-----

if(make_plots){
  
  print('making plots')
  
  #----Analysis 1: Spatial scales analyses----
  
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
  
  #----Analysis 2: Basic info during the night----
  
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
  
  #---Analysis 3: Leader/follower analyses (original version not accounting for past behavior of "other")---
  
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
  
  #Plot: leadership hierarchy (sorted by mean leading rate)
  mean_lead <- colMeans(p_approach_middle, na.rm=T)
  ord <- order(mean_lead, decreasing = F)
  plotpath <- paste0(dir, 'plots/dyadic/approach_heatmap_10-500m_sorted.png')
  png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
  image.plot(p_approach_middle[ord,ord], xaxt='n', yaxt='n' ,zlim=c(0,1),col=viridis(256),xlab='Follower',ylab='Leader', main = 'P(approach) when 10-500 m apart')
  axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las = 2)
  axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las=2)
  y <- c(matrix(rep(seq(0,1,length.out=n_inds), each = n_inds), nrow = n_inds, ncol = n_inds))
  x <- c(matrix(rep(seq(0,1,length.out=n_inds), n_inds), nrow = n_inds, ncol = n_inds))
  text(x,y,paste0(round(c(p_approach_middle[ord,ord])*100),'%'))
  dev.off()
  
  #Plots: Is the 'leadership' hierarchy consistent across nights?
  periods$date <- as.Date(periods$start_time)
  n_nights <- nrow(periods)
  p_approach_by_night <- npoints_by_night <- array(NA, dim = c(n_inds,n_inds,n_nights))
  for(p in 1:nrow(periods)){
    t0 <- periods$t0_idx[p]
    tf <- periods$tf_idx[p]
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        curr <- head_twd_night[i,j,t0:tf]
        dists_curr <- dyad_dists[i,j,t0:tf]
        idxs <- which(dists_curr >= 10 & dists_curr < 500)
        npoints_by_night[i,j,p] <- sum(!is.na(curr[idxs]))
        p_approach_by_night[i,j,p] <- mean(curr[idxs]< pi/2,na.rm=T)
      }
    }
  }
  p_approach_by_night[which(npoints_by_night < 60*10)] <- NA
  follow_by_night_mat <- apply(p_approach_by_night, c(1,3), mean, na.rm=T)
  lead_by_night_mat <- apply(p_approach_by_night, c(2,3), mean, na.rm=T)
  follow_by_night_mat[which(is.nan(follow_by_night_mat))] <- NA
  lead_by_night_mat[which(is.nan(lead_by_night_mat))] <- NA
  
  #store in a data frame
  p_approach_dat <- data.frame(follower_idx = rep(rep(1:n_inds, n_inds),n_nights), 
                               leader_idx = rep(rep(1:n_inds, each = n_inds),n_nights),
                               night_idx = rep(1:n_nights, each = n_inds*n_inds))
  p_approach_dat$prob_approach <- p_approach_by_night[cbind(p_approach_dat$follower_idx, p_approach_dat$leader_idx, p_approach_dat$night_idx)]
  p_approach_dat$follower_id <- ids$code[p_approach_dat$follower_idx]
  p_approach_dat$leader_id <- ids$code[p_approach_dat$leader_idx]
  p_approach_dat$date <- periods$date[p_approach_dat$night_idx]
  p <- ggplot(p_approach_dat, aes(x=leader_id, y=prob_approach)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2))
  
  #probability of being followed (mean across all conspecifics) each night
  lead_by_night <- data.frame(ind_idx = rep(1:n_inds, n_nights),
                              night_idx = rep(1:n_nights, each = n_inds))
  lead_by_night$lead_prob <- lead_by_night_mat[cbind(lead_by_night$ind_idx, lead_by_night$night_idx)]
  lead_by_night$follow_prob <- follow_by_night_mat[cbind(lead_by_night$ind_idx, lead_by_night$night_idx)]
  lead_by_night$code <- ids$code[lead_by_night$ind_idx]
  
  plotpath <- paste0(dir, 'plots/dyadic/lead_by_night_10-500m.png')
  p <- ggplot(lead_by_night, aes(x=code, y=lead_prob, fill = code)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_minimal() + 
    scale_fill_viridis(discrete=T) + xlab('Individual') + ylab('Mean probability of being approached')
  ggsave(plot = p, filename = plotpath)
  
  plotpath <- paste0(dir, 'plots/dyadic/follow_by_night_10-500m.png')
  p2 <- ggplot(lead_by_night, aes(x=code, y=follow_prob, fill = code)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_minimal() + 
    scale_fill_viridis(discrete=T) + xlab('Individual') + ylab('Mean probability of being approached')
  ggsave(plot = p2, filename = plotpath)
  
  #-----Analysis 4: plots for interactions analysis (accounting for past behavior of "other")----
  
  #Separate out probability of approach when in different situations:
  #1. "Other" is moving away from the ego - 'follow'
  #2. "Other" is moving toward from the ego - 'converge'
  #3. "Other" is stationary - 'join'
  #Two versions - either we include the ego being stationary as a non-response (_all), or we exclude stationary cases and therefore condition on the ego moving (_moving)
  
  #Aggregate patterns across dyads
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/toward_all.png'), other_behav = NULL, include_moving_only = F, stratify_by = 'ego')
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/toward_moving.png'), other_behav = NULL, include_moving_only = T)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_all.png'), other_behav = 'away', include_moving_only = F)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_moving.png'), other_behav = 'away', include_moving_only = T)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/converge_all.png'), other_behav = 'toward', include_moving_only = F)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/converge_moving.png'), other_behav = 'toward', include_moving_only = T)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/join_all.png'), other_behav = 'stationary', include_moving_only = F)
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/join_moving.png'), other_behav = 'stationary', include_moving_only = T)

  #Comparing patterns across different dyads
  
  #Following
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_moving_by_ego.png'), other_behav = 'away', include_moving_only = T, stratify_by = 'ego')
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_moving_by_other.png'), other_behav = 'away', include_moving_only = T, stratify_by = 'other')
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_all_by_ego.png'), other_behav = 'away', include_moving_only = F, stratify_by = 'ego')
  make_approach_vs_dist_plot(approach_data, dist_bins, plotpath = paste0(dir, 'plots/interaction_types/follow_all_by_other.png'), other_behav = 'away', include_moving_only = F, stratify_by = 'other')
  
  #Could add in analyses of converging and joining as well, but will skip it for now
  
  #"Following" heat map for the most relevant distance bin (5 - 200 m)
  p_approach_middle <- n_middle <- matrix(NA, nrow = n_inds, ncol = n_inds)
  dat_curr <- approach_data[which(approach_data$past_behav_j=='away'),]
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      idxs <- which(dat_curr$dyad_dist >= 5 & dat_curr$dyad_dist < 200 & dat_curr$i_ego == i & dat_curr$j_other == j)
      p_approach_middle[i,j] <- mean(dat_curr$i_approaches_j_given_moving[idxs], na.rm=T)
      n_middle[i,j] <- sum(!is.na(dat_curr$i_approaches_j_given_moving[idxs]))
    }
  }
  
  mean_lead <- colMeans(p_approach_middle, na.rm=T)
  ord <- order(mean_lead, decreasing = F)
  
  plotpath <- paste0(dir, 'plots/interaction_types/follow_moving_heatmap_5-200m.png')
  png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
  image.plot(p_approach_middle[ord,ord], xaxt='n', yaxt='n' ,zlim=c(min(p_approach_middle,na.rm=T),1),col=viridis(256),xlab='Follower',ylab='Leader', main = 'P(follow) when 5-200 m apart')
  axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las = 2)
  axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las=2)
  y <- c(matrix(rep(seq(0,1,length.out=n_inds), each = n_inds), nrow = n_inds, ncol = n_inds))
  x <- c(matrix(rep(seq(0,1,length.out=n_inds), n_inds), nrow = n_inds, ncol = n_inds))
  text(x,y,paste0(round(c(p_approach_middle[ord,ord])*100),'%'))
  dev.off()
  
  #Could add heatmaps for other behaviors, but skipping for now
  
}



#TODO:
#look into autocorrelation, consider downsampling 
#model with autocorrelation model to get error bars 
#nearest neighbor
#different time scales - 5 sec, 10 sec, 30, sec, 1 min ,5, 10

