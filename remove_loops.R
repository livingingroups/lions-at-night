#exploring how to remove loops that are likely gps errors

#For each time point, we check whether the individual moved at least a distance R_outer away during the subsequent
#max_time_steps time steps. If so, we then check how long it was until that individual returned to the original
#location (i.e. within a distance R_inner of that location). This identifies a loop.

#We then compute various metrics describing these loops, including:
#loops$dt: duration in sec
#loops$dist_traveled: total distance traveled
#loops$max_dist: maximum distance away from the original point 
#loops$speed_mean: mean speed (m/s)
#loops$vedba_mean: mean vedba value during the loop
#loops$vedba_max: max vedba during the loop
#loops$vedba_median: median vedba during the loop
#loops$vedba_95pct: 95% quantile of vedba during the loop
#loops$active_pct: percent of time the hyena was "active" (above vedba_thresh)

library(rhdf5)
library(lubridate)

options(digits=6)

#PARAMS - adjustable

#Parameters for identifying loops
R_inner <- 20 #maximum distance from the first point to complete a loop
R_outer <- 300 #minimum distance away from the first point to count as a loop
max_time_steps <- 6 #maximum number of time steps it takes to go from the initial locaiton to a distance R_outer away, to start a loop

#Parameters for finding short, fast, low vedba loops
max_loop_duration <- 20*60 #maximum number of time steps for an identified loop to be considered "short"
vedba_thresh <- -3.4 #threshold value of vedba to be considered in "low activity" state (-3.4 from Minasandra et al. 2023)
min_loop_log_speed <- -4 #threshold value of median

#PARAMS - keep fixed
sec_per_timestep <- 30 #seconds per timestep in gps data (30)
acc_dir <- '~/EAS_shared/hyena/working/processed/cc23/Acc_Data/'
acc_fs <- 25 #sample rate of ACC data (25 Hz)


#----LOAD DATA

#--GPS data
load('~/EAS_shared/hyena/working/processed/cc23/hyena_south2023_xy_level1_30s.RData')

#number of individuals and time steps
n_inds <- nrow(xs)
n_times <- ncol(xs)

#--ACC data
#get acc data and store in a matrix vedbas
setwd(acc_dir)
files <- list.files()
n_times_acc <- as.numeric(difftime(timestamps[n_times],timestamps[1],units='secs')) * acc_fs
vedba_timestamps <- seq.POSIXt(from = timestamps[1], to = timestamps[n_times], length.out = n_times_acc)
vedbas <- matrix(NA, nrow = n_inds, ncol = n_times_acc)
first_time <- as.numeric(vedba_timestamps[1])
print('getting vedba data')
for(i in 1:n_inds){
  print(i)
  id <- ids$id_code[i]
  filename <- paste0('cc23_',toupper(id),'_A_25Hz.h5')
  if(filename %in% files){
    f <- rhdf5::H5Fopen(filename)
    vedba <- rhdf5::h5read(f, 'vedba')
    utc_vec <- rhdf5::h5read(f, 'UTC')
    utc <- as.POSIXct(strftime(paste0(utc_vec[1],'-',utc_vec[2],'-',utc_vec[3],' ',utc_vec[4],':',utc_vec[5],':',utc_vec[6]),format='%Y-%m-%d %H:%M:%OS4'),tz='UTC')
    first_time_idx <- (as.numeric(utc) - first_time) * acc_fs + 1
    
    vedbas[i, first_time_idx:(first_time_idx + length(vedba)-1)] <- vedba
  }
}

#----FIND LOOPS

print('finding loops')
loops <- data.frame()

#loop over individuals
for(ind in 1:n_inds){
  xi <- xs[ind,]
  yi <- ys[ind,]
  
  #for each time step, check if it starts a loop
  for(t in 1:(n_times-max_time_steps)){
    
    x0 <- xi[t]
    y0 <- yi[t]
    
    #skip if NA
    if(is.na(x0)){
      next
    }
    
    #get the time the inner radius was crossed, if within the max_time_steps time threshold
    t_R <- NA
    for(t_next in (t+1):(t+max_time_steps)){
      dx <- xi[t_next] - x0
      dy <- yi[t_next] - y0
      dist <- sqrt(dx^2 + dy^2)
      if(is.na(dist)){
        break
      }
      if(dist >= R_outer){
        t_R <- t_next
        break
      }
    }
    
    #find the next subsequent time the dist is within R_inner of the original location (t_return)
    if(!is.na(t_R)){
      
      t_return <- NA
      for(t_next in (t_R+1):n_times){
        dx <- xi[t_next] - x0
        dy <- yi[t_next] - y0
        dist <- sqrt(dx^2 + dy^2)
        if(is.na(dist)){
          break
        }
        if(dist < R_inner){
          t_return <- t_next
          loops <- rbind(loops, data.frame(ind = ind, t0 = t, t_R = t_R, t_return = t_return))
          break
        }
      }
    }
  }
}

#remove events that are contiguous with other events
dt <- diff(loops$t0)
loops <- loops[which(dt>1),]


print('getting loop characteristics')
#get time until return (duration of the loop)
loops$dt <- (loops$t_return - loops$t0) * sec_per_timestep

#initial time in UTC
loops$t0_UTC <- timestamps[loops$t0]

#total distance traveled during the loop, and max distance away from the initial location
loops$dist_traveled <- loops$max_dist <- NA
for(i in 1:nrow(loops)){
  t0 <- loops$t0[i]
  tf <- loops$t_return[i]
  ind <- loops$ind[i]
  
  xi <- xs[ind, t0:tf]
  yi <- ys[ind, t0:tf]
  
  dists <- sqrt( (xi - xi[1])^2 + (yi - yi[1])^2  )
  disps <- sqrt(diff(xi)^2 + diff(yi)^2)
  
  loops$dist_traveled[i] <- sum(disps, na.rm=T)
  loops$max_dist[i] <- max(dists, na.rm=T)
  
  
}

#get average speed
loops$speed_mean <- loops$dist_traveled / (loops$dt*sec_per_timestep)

#get mean, median, max, and 95 percentile of vedba during the loop
loops$t0_idx_vedba <- loops$t_return_idx_vedba <- NA
loops$vedba_mean <- loops$vedba_median <- loops$vedba_max <- loops$vedba_95pct <- loops$pct_active <- NA
for(i in 1:nrow(loops)){
  
  #get time indexes into vedba matrix
  ind <- loops$ind[i]
  loops$t0_idx_vedba[i] <- (loops$t0[i]-1)*acc_fs*sec_per_timestep + 1
  loops$t_return_idx_vedba[i] <- (loops$t_return[i]-1)*acc_fs*sec_per_timestep + 1
  
  #get mean and max vedba
  vedba_curr <- vedbas[ind, loops$t0_idx_vedba[i]:loops$t_return_idx_vedba[i]]
  loops$vedba_mean[i] <- mean(vedba_curr, na.rm=T)
  loops$vedba_median[i] <- median(vedba_curr, na.rm=T)
  loops$vedba_max[i] <- max(vedba_curr, na.rm=T)
  loops$vedba_95pct[i] <- quantile(vedba_curr, 0.95, na.rm=T)
  loops$pct_active[i] <- mean(log(vedba_curr) > vedba_thresh, na.rm=T)
  
}

#individual id
loops$id <- ids$id_code[loops$ind]

#indexes to short duration loops, low vedba loops, high speed loops
short <- which(loops$dt < max_loop_duration)
low_vedba <- which(loops$pct_active < .6)
high_speed <- which(log(loops$speed_mean) > min_loop_log_speed)

#short duration and high speed loops
short_and_fast <- intersect(high_speed, short)

#short duration and high speed with low percent active
short_and_fast_and_inactive <- intersect(short_and_fast, low_vedba)

#----PLOTTING EXAMPLES

plot_loop_event <- function(xs, ys, loops, i){
  par(mfrow=c(2,1),mar=c(2,0,0,0))
  plot(xs[loops$ind[i],(loops$t0[i]-60):(loops$t_return[i]+60)],ys[loops$ind[i],(loops$t0[i]-60):(loops$t_return[i]+60)],asp=1,type='l',col='black',xlab='',ylab='')
  lines(xs[loops$ind[i],(loops$t0[i]-60):loops$t0[i]],ys[loops$ind[i],(loops$t0[i]-60):loops$t0[i]],asp=1,type='l',col='blue')
  lines(xs[loops$ind[i],loops$t0[i]:(loops$t_return[i])],ys[loops$ind[i],loops$t0[i]:(loops$t_return[i])],col='red')
  lines(xs[loops$ind[i],loops$t_return[i]:(loops$t_return[i]+60)],ys[loops$ind[i],loops$t_return[i]:(loops$t_return[i]+60)],col='darkgreen')
  points(xs[loops$ind[i],loops$t0[i]:loops$t_return[i]],ys[loops$ind[i],loops$t0[i]:loops$t_return[i]],asp=1,col='red',cex=.4)
  vedba_idxs <- loops$t0_idx_vedba[i]:(loops$t0_idx_vedba[i]+loops$dt[i]*acc_fs)
  plot(seq(0,(length(vedba_idxs)-1)/acc_fs, 1/acc_fs)/60,log(vedbas[loops$ind[i],vedba_idxs]),type='l',xlab='Samples',ylab='Vedba',ylim=c(-8,2))
  abline(h=vedba_thresh,col='red')
}

plot_loop_event(xs, ys, loops, short_and_fast_and_inactive[1])

    

