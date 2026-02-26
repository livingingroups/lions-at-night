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
dist_bins <- c(0, 10^seq(0.6,3,.1)) #distance between dyads bins to use for probability of approaching analyses
min_points_to_plot <- 3*60*60 #minimum of 3 hrs of data needed to plot things (currently only used in analysis 3 and 4)

#parameters for spatial scales analyses
#this uses windows centered on a given value and going val - pct_win/100*val to val + pct_win/100*val
min_dist_bin <- 5 #minimum distance bin
max_dist_bin <- 2000 #maximum distance bin
n_bins <- 50 #number of bins (log spaced)
pct_win <- 40 #percentage window for each bin (bin will run for x - pct_win*x/100 to x + pct_win*x/100)
n_boots <- 100 #number of bootstraps to do (resampling nights) to produce error bars
downsamp_rate <- 10 #downsampling rate for approach data (to save computation)

calculate_approach_data <- T #whether the generate the data frame of approach data
load_approach_data <- F #whether the load the cohesion dataframe from file
savename_df <- 'data_namibia/processed/cohesion_data.RData' #where to save the data frame 'approach_data' (from analysis 4) - if NULL, nothing will be saved

make_cohesion_plots <-  F #whether to generate cohesion (spatial scales) plots or not
make_leadership_plots <- F #whether to generate the leadership / influence plots or not

#-------------FUNCTIONS----------------

#helper function to get relevant statistics of the distribution of values for a set of bins (binned by another variable x)
#vals = values for which you are getting the distribution of (e.g. heading correlation)
#x = what you are binning on (e.g. dyadic distance)
#bins = bins to use for binning x
#subset_ids = identity of the non-independent subset (e.g. day) associated with each value, to use for bootstrapping
#the dimensions of vals, x, and subset_ids must be the same
get_distrib_statistics <- function(vals, x, subset_ids = NULL, n_boots = 100, min_dist=5, max_dist=1500, n_bins=100, pct_win=50){
  out <- list()
  
  #logarithmically spaced bin midpoints
  bins <- exp(seq(log(min_dist), log(max_dist), length.out = n_bins))
  out$bins <- bins
  out$mean <- rep(NA, length(bins)-1)
  #loop over all bins
  for(i in 1:length(bins)){
    
    #get indexes to data in that bin
    idxs <- which(x >= bins[i]*(1-pct_win/100) & x < bins[i]*(1+pct_win/100))
    
    #get relevant stats
    out$mean[i] <- mean(vals[idxs], na.rm=T)
    
    out$n[i] <- sum(!is.na(vals[idxs]))
  }
  
  if(!is.null(subset_ids)){
    subsets <- unique(c(subset_ids))
    subsets <- subsets <- subsets[which(!is.na(subsets))]
    n_subsets <- length(subsets)
    means_boot <- matrix(NA, nrow = n_boots, ncol = n_bins)
    for(b in 1:n_boots){
      #construct a bootstrapped dataset by drawing with replacement from the subsets
      subsets_to_use <- sample(subsets, replace = T)
      x_boot <- vals_boot <- c()
      for(s in subsets_to_use){
        subset_idxs <- which(subset_ids == s)
        x_boot <- c(x_boot, x[subset_idxs])
        vals_boot <- c(vals_boot, vals[subset_idxs])
      }
      
      #measure the mean of the bootstrapped dataset by bin
      for(i in 1:length(bins)){
        
        #get indexes to data in that bin
        idxs <- which(x_boot >= bins[i]*(1-pct_win/100) & x_boot < bins[i]*(1+pct_win/100))
        means_boot[b,i] <- mean(vals_boot[idxs], na.rm=T)
      }
      
    }
    out$means_boot <- means_boot
  }
  
  return(out)
}

#helper function to make plots of different metrics vs distance
#plot include mean (thick black lines) and bootstrapped replicates (based on subset_ids, if present)
#INPUTS:
# dat: list containing 
#   dat$bins: numeric vector of bin centers
#   dat$mean: numeric vector of values to plot (thick lines)
#   dat$means_boot: [n_boots x n_bins] dimension matrix of bootstrapped replicate data
make_plot <- function(dat, plotpath, ylab, ylim = NULL, xlab = 'Distance apart (m)', abline_y = NULL, logx = F){
  
  if(is.null(ylim)){
    ylim <- c(min(dat$means_boot,na.rm=T), max(dat$means_boot, na.rm=T))
  }
  
  png(filename = plotpath, width = 6, height = 6, units = 'in', res = 300)
  mids <- dat$bins
  if(logx){
    plot(mids, dat$mean, xlab = xlab, ylab = ylab, pch = 19, col = '#00000000', cex = 2, ylim = ylim ,log='x', cex.axis = 1.5, cex.lab = 1.5)
  } else{
    plot(mids, dat$mean, xlab = xlab, ylab = ylab, pch = 19, col = '#00000000', cex = 2, ylim = ylim, cex.axis = 1.5, cex.lab = 1.5)
    
  }
  if(!is.null(abline_y)){
    abline(h=abline_y, lty = 2, lwd = 2)
  }
  
  #plot bootstrap replicates
  for(b in 1:nrow(dat$means_boot)){
    lines(mids, dat$means_boot[b,], col = '#00000011', lwd = 1)
  }
  #plot mean
  lines(mids, dat$mean, col = 'black', lwd = 3)

  dev.off()
}

#Plotting function for 'interactions' (i.e. following, converging, joining)
#other_behav should be either NULL (use all data), 'toward', 'away', or 'stationary'
#stratify_by gives whether to stratify by dyad, can be 'ego','other' or NULL (if you want to make the plot aggregating across all dyads)
make_approach_vs_dist_plot <- function(approach_data, plotpath, other_behav = NULL, include_moving_only = F, stratify_by = NULL, min_points_to_plot = 60*60,
                                       min_dist_bin = 5, max_dist_bin = 1000, n_bins = 50, pct_win = 40){
  
  #if indexes to subset by aren't specified, use all rows of data
  if(is.null(other_behav)){
    idxs_to_use <- 1:nrow(approach_data)
  } else{
    idxs_to_use <- which(approach_data$past_behav_j == other_behav)
  }
  
  #subset to only data we need to use
  dat <- approach_data[idxs_to_use,]
  
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
      toward <- get_distrib_statistics(vals = dat$i_approaches_j_given_moving, x = dat$dyad_dist, 
                                       min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
    } else{
      toward <- get_distrib_statistics(vals = dat$i_approaches_j, x = dat$dyad_dist, 
                                       min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
    }
    
    #make plot
    png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
    if(include_moving_only){
      plot(toward$bins, toward$mean, pch = 19, cex = 2, type = 'l', lwd = 3, log = 'x', xlab = 'Distance apart (m)', ylab = 'Probability of approach (given moved)', main = main_lab)
      abline(h=0.5, lty = 2)
    } else{
      plot(toward$bins, toward$mean, pch = 19, cex = 2, type = 'l', lwd = 3, log = 'x', xlab = 'Distance apart (m)', ylab = 'Probability of approach', main = main_lab)
    }
    dev.off()
  } else{
    n_inds <- max(approach_data$i_ego)
    approach_probs <- ns <- array(NA, dim = c(n_inds, n_inds, n_bins))
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        if(i!=j){
          idxs <- which(dat$i_ego==i & dat$j_other==j)
          #get means by distance bin
          if(include_moving_only){
            toward <- get_distrib_statistics(vals = dat$i_approaches_j_given_moving[idxs], x = dat$dyad_dist[idxs], 
                                             min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
          } else{
            toward <- get_distrib_statistics(vals = dat$i_approaches_j[idxs], x = dat$dyad_dist[idxs], 
                                             min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
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
        plot(NULL, xlim = c(min(toward$bins), max(toward$bins)), ylim = ylims, log = 'x', main = paste0('Leader = ',ids$code[j]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
        for(i in 1:n_inds){
          if(i!=j){
            curr <- approach_probs[i,j,]
            lines(toward$bins, curr, lwd = 2, col = cols[i])
          }
        }
        if(include_moving_only){
          abline(h=0.5, lty = 2)
        }
        legend('bottomleft', col = cols, lwd = 2, legend = paste0('Follower = ',ids$code))
      }
    }
    if(stratify_by == 'ego'){
      for(i in 1:n_inds){
        plot(NULL, xlim = c(min(toward$bins), max(toward$bins)), ylim = ylims, log = 'x', main = paste0('Follower = ',ids$code[i]), ylab = 'P(approach)',xlab = 'Distance apart (m)')
        for(j in 1:n_inds){
          if(i!=j){
            curr <- approach_probs[i,j,]
            lines(toward$bins, curr, lwd = 2, col = cols[j])
          }
        }
        if(include_moving_only){
          abline(h=0.5, lty = 2)
        }
        legend('bottomleft', col = cols, lwd = 2, legend = paste0('Leader = ',ids$code))
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

#Get periods of data to use for spatial scales and leadership analyses (only nighttime data)
start_times <- as.POSIXct(paste(seq.Date(date('2023-05-06'),date('2023-06-15'),by=1), '16:00:00'), tz = 'UTC') #timestamps of times to start each analysis period (night) - replace eventually with info on when the group left its rest spot
end_times <- as.POSIXct(paste(seq.Date(date('2023-05-07'),date('2023-06-16'),by=1), '06:00:00'), tz = 'UTC') #timestamps of times to end each analysis period (night) - replace eventually with info on when the group left its rest spot

periods <- data.frame(start_time = start_times, end_time = end_times)
periods$t0_idx <- match(periods$start_time, timestamps)
periods$tf_idx <- match(periods$end_time, timestamps)
periods$subset_id <- seq(1, nrow(periods))

#get sequence of subset ids (e.g. nights or periods)
#this identifies each data point with a period that it is part of, matched to subset_id column of periods data frame
subset_ids <- array(NA, dim = c(n_inds,n_inds,n_times))
for(i in 1:nrow(periods)){
  subset_ids[,,periods$t0_idx[i]:periods$tf_idx[i]] <- periods$subset_id[i]
}

#relevant time indexes for spatial scales analyses
relevant_t_idxs <- c()
for(i in 1:nrow(periods)){
  relevant_t_idxs <- c(relevant_t_idxs, periods$t0_idx[i]:periods$tf_idx[i])
}
idxs_to_remove <- setdiff(1:n_times, relevant_t_idxs)


#COMPUTE METRICS

#headings and speeds (future)
print('...headings and speeds...')
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
print('...dyadic distances...')
dyad_dists <- cocomo::get_group_dyadic_distances(xs, ys)

#group elongation
print('...elongation...')
out <- cocomo::get_group_elongation(xs, ys, min_inds_tracked = n_inds - 1)
elongation <- out$elongation
long_axis_angle <- out$long_axis_angle

#group polarization (computed for all instances where at least 5 individuals were tracked, otherwise NA)
print('...polarization...')
polarization <- cocomo::get_group_polarization(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)

#group speed (computed for all instances where at least 5 individuals were tracked, otherwise NA)
print('...group speed...')
out <- cocomo::get_group_heading_and_speed(xs, ys, heading_type = 'temporal', t_window = speed_dt, min_inds_tracked = n_inds-1)
group_speed <- out$speeds
group_head <- out$heads


#get mean dyadic distance (computed for all instances where at least 5 individuals were tracked, otherwise NA)
mean_dyad_dist <- apply(dyad_dists, 3, FUN = mean, na.rm=T)
n_tracked <- colSums(!is.na(xs))
mean_dyad_dist[which(n_tracked < n_inds - 1)] <- NA

#dyadic distance changes, speed differences, and heading correlations between pairs of individuals
print('...dyadic distance diffs, speed diffs, heading corrs...')
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
print('...relative headings...')
head_twd <- array(NA, dim = c(n_inds, n_inds, n_times))
for(i in 1:n_inds){
  for(j in 1:n_inds){
    if(i != j){
      
      #heading of i
      head_i <- heads[i,] 
      
      #vector pointing from i to j - check this bit Gen
      dx_ij <- (xs[j,] - xs[i,]) / dyad_dists[i,j,]
      dy_ij <- (ys[j,] - ys[i,]) / dyad_dists[i,j,]
      
      #angle between these vectors - check this bit Gen
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
#check this gen
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

if(calculate_approach_data){
  print('generating cohesion dataframe (approach_data)')
  #----Data frame to hold all 'interactions' data at night----
  #i is the "ego", i.e. the one whose behavior in response to "other" we are considering - follower
  #j is the "other, i.e. the reference individual relative to whom i's behavior is being considered - leader
  approach_data <- data.frame(i_ego = rep(rep(1:n_inds, n_inds),n_times), 
                              j_other = rep(rep(1:n_inds, each = n_inds),n_times),
                              tidx = rep(1:n_times, each = n_inds*n_inds))
  
  #remove non-night columns (otherwise we might run out of memory :/)
  approach_data <- approach_data[which(approach_data$tidx %in% relevant_t_idxs),]
  
  #downsample to only indexes that are divisible by downsamp_rate = 10 (save computational power)
  approach_data <- approach_data[which(approach_data$tidx %% downsamp_rate == 0),]
  
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
  #elongation - elongation of the group at time step tidx
  #long_axis_angle - angle of the long axis of the group at time step tidx
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
  
  #group metrics - polarization, elongation, long axis angle, group heading, group speed, group tilt
  approach_data$group_polarization <- polarization[approach_data$tidx]
  approach_data$group_elongation <- elongation[approach_data$tidx]
  approach_data$group_long_axis_angle <- long_axis_angle[approach_data$tidx]
  approach_data$group_heading <- group_head[approach_data$tidx]
  approach_data$group_speed <- group_speed[approach_data$tidx]
  approach_data$group_mean_dyad_dist <- mean_dyad_dist[approach_data$tidx]
  tilt_ang <- acos(cos(approach_data$group_long_axis_angle)*cos(approach_data$group_heading) + 
                                     sin(approach_data$group_long_axis_angle)*sin(approach_data$group_heading))
  tilt_ang[which(tilt_ang > pi/2)] <- pi - tilt_ang[which(tilt_ang > pi/2)]
  approach_data$group_tilt <- tilt_ang / (pi/2)
  
  #important bit - gen check
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
} else{
  if(load_approach_data){
    print('loading cohesion data')
    load(file=paste0(dir, savename_df))
  }
}

#-----PLOTS-----

if(make_cohesion_plots){
  
  print('making cohesion plots')
  
  #----Analysis 1: Spatial scales analyses----
  
  #calculate mean, median, iqr of metrics during relevant periods only, as a function of distance apart
  metrics <- list()
  metrics$dist_bins <- dist_bins
  
  #get statistics of the distribution for each metric within each bin
  downsamp_idxs <- seq(1, n_times, downsamp_rate) #downsample to save computational time
  metrics$speed_diff <- get_distrib_statistics(speed_diffs_night[,,downsamp_idxs], dyad_dists[,,downsamp_idxs], subset_ids[,,downsamp_idxs], n_boots = n_boots, min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
  metrics$dyad_dist_change <- get_distrib_statistics(dyad_dist_changes_night[,,downsamp_idxs], dyad_dists[,,downsamp_idxs], subset_ids[,,downsamp_idxs], n_boots = n_boots, min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
  metrics$ang_between_heads <- get_distrib_statistics(acos(head_corrs_night[,,downsamp_idxs])*180/pi, dyad_dists[,,downsamp_idxs], subset_ids[,,downsamp_idxs], n_boots = n_boots, min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
  metrics$prob_approach <- get_distrib_statistics(head_twd_night[,,downsamp_idxs]<=pi/2, dyad_dists[,,downsamp_idxs], subset_ids[,,downsamp_idxs], n_boots = n_boots, min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)

  #Plot: Angle between headings 
  plotpath <- paste0(dir, 'plots/spatial_scales/ang_between_heads.png')
  ylab <- "Mean angle between headings (degrees)"
  ylim <- c(0,120)
  dat <- metrics$ang_between_heads
  make_plot(dat, plotpath, ylab, ylim, abline_y = 90, logx = T)
  
  #Plot: Mean difference in speed as a function of distance 
  plotpath <- paste0(dir, 'plots/spatial_scales/abs_speed_diff.png')
  ylab = 'Mean absolute speed difference (m/s)'
  dat <- metrics$speed_diff
  make_plot(dat, plotpath, ylab, logx = T, ylim = c(.15,.7))
  
  #Plot: Change in dyadic distance
  plotpath <- paste0(dir, 'plots/spatial_scales/dyad_dist_change.png')
  ylab = 'Mean change in dyadic distance (m)'
  dat <- metrics$dyad_dist_change
  make_plot(dat, plotpath, ylab, abline_y = 0, logx = T, ylim = c(-1,1))
  
  #Plot: Change in dyadic distance
  plotpath <- paste0(dir, 'plots/spatial_scales/p_approach.png')
  ylab = 'Probability of approach'
  dat <- metrics$prob_approach
  make_plot(dat, plotpath, ylab, abline_y = 0.5, logx = T, ylim = c(0.4,0.6))
  
  #----Analysis 2: Basic info during the night----
  
  #Plot: Distribution of log(dyadic distances) between individuals
  png(filename = paste0(dir, 'plots/overall_hists/hist_log_dyad_dist.png'), width = 8, height = 6, units = 'in', res = 300)
  histo <- hist(log(dyad_dists, 10), plot=T, breaks=80, xlab = 'Dyadic distance (m) - log bins',main='', xlim = c(0,4), freq = F, xaxt ='n')
  axis(1, at = seq(0,4,1), labels = 10^seq(0,4,1))
  dev.off()
  
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
  
}

if(make_leadership_plots){
  
  #-----Dyad-level leadership plots----
  
  #For the purposes of this analysis, we define the probability of following as:
  #Probability of an individual i (the 'follower') moving toward an individual j (the 'leader') given the individual j has, in the past, moved away from i
  #Numerator: past_behav_j == 'away' && fut_behav_i == 'toward'
  #Denominator = past_behav_j == 'away' & (fut_behav_i == 'toward' or fut_behav_i = 'away')
  #In this case, if the follower remains stationary, it is not counted as not following (i.e. not included in the denominator)
  #Another way to look at it is we are only analyzing instances in which both individuals are moving
  #Conditioned on that, if i moves away from j, what does j do (does it move toward i or move away from i?)
  
  #Plot 1: Do this across all dyads aggregated across all data
  #make_approach_vs_dist_plot(approach_data_downsamp, plotpath = paste0(dir, 'plots/leader_follower/p_follow_all_inds.png'), other_behav = 'away', include_moving_only = F, stratify_by = NULL, min_points_to_plot = 60)
  make_approach_vs_dist_plot(approach_data, plotpath = paste0(dir, 'plots/leader_follower/p_follow_all_inds_moving.png'), other_behav = 'away', include_moving_only = T, stratify_by = NULL, min_points_to_plot = 60)
  
  #Plot 2: Stratify by 'other' (leader)
  #make_approach_vs_dist_plot(approach_data_downsamp, plotpath = paste0(dir, 'plots/leader_follower/p_follow_by_leader.png'), other_behav = 'away', include_moving_only = F, stratify_by = 'other', min_points_to_plot = 60)
  make_approach_vs_dist_plot(approach_data, plotpath = paste0(dir, 'plots/leader_follower/p_follow_by_leader_moving.png'), other_behav = 'away', include_moving_only = T, stratify_by = 'other', min_points_to_plot = 60)
  
  #Plot 3: Stratify by 'ego' (follower)
  make_approach_vs_dist_plot(approach_data, plotpath = paste0(dir, 'plots/leader_follower/p_follow_by_follower_moving.png'), other_behav = 'away', include_moving_only = T, stratify_by = 'ego', min_points_to_plot = 60)
  
  #Plot 4: "Following" heat map for the most relevant distance bin (0 - 200 m)
  # p_approach_middle <- n_middle <- matrix(NA, nrow = n_inds, ncol = n_inds)
  # dat_curr <- approach_data[which(approach_data$past_behav_j=='away'),]
  # for(i in 1:n_inds){
  #   for(j in 1:n_inds){
  #     idxs <- which(dat_curr$dyad_dist >= 0 & dat_curr$dyad_dist < 300 & dat_curr$i_ego == i & dat_curr$j_other == j)
  #     p_approach_middle[i,j] <- mean(dat_curr$i_approaches_j_given_moving[idxs], na.rm=T)
  #     n_middle[i,j] <- sum(!is.na(dat_curr$i_approaches_j_given_moving[idxs]))
  #   }
  # }
  # 
  # mean_lead <- rowMeans(p_approach_middle, na.rm=T)
  # ord <- order(mean_lead, decreasing = T)
  # 
  # plotpath <- paste0(dir, 'plots/leader_follower/follow_moving_heatmap_5-200m.png')
  # png(filename = plotpath, width = 8, height = 6, units = 'in', res = 300)
  # image.plot(p_approach_middle[ord,ord], xaxt='n', yaxt='n' ,zlim=c(min(p_approach_middle,na.rm=T),1),col=viridis(256),xlab='Follower',ylab='Leader', main = 'Probability of following when moving (<300 m apart)')
  # axis(1, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las = 1)
  # axis(2, at = seq(0,1,length.out=n_inds), labels = ids$code[ord], las=2)
  # y <- c(matrix(rep(seq(0,1,length.out=n_inds), each = n_inds), nrow = n_inds, ncol = n_inds))
  # x <- c(matrix(rep(seq(0,1,length.out=n_inds), n_inds), nrow = n_inds, ncol = n_inds))
  # text(x,y,paste0(round(c(p_approach_middle[ord,ord])*100),'%'))
  # dev.off()
  
  plotpath <- paste0(dir, 'plots/leader_follower/follow_moving_heatmap_5-300m.png')
  make_combined_heatmap_lines_plot(approach_data = approach_data, plotpath = plotpath, other_behav = 'away',include_moving_only = T, min_points_to_plot = 100, 
                                   min_dist_bin = 5, max_dist_bin = 300, n_bins = 100, pct_win = 40)
  
}


make_combined_heatmap_lines_plot <- function(approach_data, plotpath, other_behav = NULL, include_moving_only = F, min_points_to_plot = 100,
                                       min_dist_bin = 5, max_dist_bin = 300, n_bins = 50, pct_win = 40){
  
  #if indexes to subset by aren't specified, use all rows of data
  if(is.null(other_behav)){
    idxs_to_use <- 1:nrow(approach_data)
  } else{
    idxs_to_use <- which(approach_data$past_behav_j == other_behav)
  }
  
  #subset to only data we need to use
  dat <- approach_data[idxs_to_use,]
  
  n_inds <- max(approach_data$i_ego)
  approach_probs <- ns <- array(NA, dim = c(n_inds, n_inds, n_bins))
  mean_probs <- matrix(NA, nrow = n_inds, ncol = n_inds)
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      if(i!=j){
        idxs <- which(dat$i_ego==i & dat$j_other==j)
        #get means by distance bin
        if(include_moving_only){
          toward <- get_distrib_statistics(vals = dat$i_approaches_j_given_moving[idxs], x = dat$dyad_dist[idxs], 
                                           min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
        } else{
          toward <- get_distrib_statistics(vals = dat$i_approaches_j[idxs], x = dat$dyad_dist[idxs], 
                                           min_dist = min_dist_bin, max_dist = max_dist_bin, n_bins = n_bins, pct_win = pct_win)
        }
        approach_probs[i,j,] <- toward$mean
        ns[i,j,] <- toward$n
        mean_probs[i,j] <- mean(dat$i_approaches_j_given_moving[idxs], na.rm=T)
      }
    }
  }
  approach_probs[which(ns < min_points_to_plot)] <- NA
  bins <- toward$bins
  ylims <- c(0.5, 1)
  xlims <- c(min(bins),max(bins))
  col_frac_matrix <- (mean_probs - min(mean_probs,na.rm=T)) / (max(mean_probs,na.rm=T) - min(mean_probs,na.rm=T))
  cols <- viridis(1001)
  
  ord <- order(rowMeans(mean_probs, na.rm=T), decreasing = F)
  
  #get center of plot
  x_center <- exp(mean(log(xlims)))
  y_center <- mean(ylims)
  
  #make the plot
  png(filename = plotpath, width = 8, height = 8, units = 'in', res = 300)
  par(mfcol = c(n_inds,n_inds), mar = c(0,0,0,0))
  for(i in ord){
    for(j in ord){
      
      #get color
      col_ij <- cols[round(col_frac_matrix[i,j]*1000)+1]
      col_ij <- paste0(substr(col_ij, 1, 7), 'CC')
    
      plot(NULL, xlim = c(min(bins), max(bins)), ylim = ylims, yaxt = 'n',xaxt='n',ylab='',xlab='',log='x', type = 'l', bty = 'n')
      print(xlims)
      print(ylims)
      if(i!=j){
        polygon(x = c(xlims[1],xlims[2],xlims[2],xlims[1]),y = c(ylims[1],ylims[1],ylims[2],ylims[2]), col = col_ij, border=NA)
        lines(bins, approach_probs[i,j,], lwd = 3)
        pct_text <- paste0(round(mean_probs[i,j]*100),'%')
        text(x_center, y_center, pct_text, cex = 1.5)
        text(xlims[1],ylims[1], min(bins))
        text(xlims[2],ylims[1], max(bins))
        text(xlims[1],ylims[2], labels = ylims[2])
        text(xlims[1], ylims[1], labels = ylims[1])
        lines(c(xlims[1], xlims[1]), c(ylims[1], ylims[2]))
        lines()
      } else{
        polygon(x = c(xlims[1],xlims[2],xlims[2],xlims[1]),y = c(ylims[1],ylims[1],ylims[2],ylims[2]), col = 'gray', border=NA)
        text(x_center, y_center, ids$code[i], cex = 2)
        
      }
    }
  }
  dev.off()
}



