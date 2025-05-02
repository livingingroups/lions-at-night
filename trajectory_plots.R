#Make some plots of trajectories

load('~/Dropbox/lions_at_night/data_namibia/processed/lion_omukutu_xy_level1.RData')


#times to plot
day <- 2
step <- 60

t0 <- 24*60*60*day - 12*60*60 + 1
tf <- t0 + 24*60*60


n_inds <- nrow(xs)
ind_cols <- c('red','orange','darkgreen','blue','purple','black')

xmin <- min(xs[,t0:tf], na.rm=T)
xmax <- max(xs[,t0:tf], na.rm=T)
ymin <- min(ys[,t0:tf], na.rm=T)
ymax <- max(ys[,t0:tf], na.rm=T)

plot(NULL, xlim=c(xmin, xmax), ylim = c(ymin, ymax),asp=1)
for(i in 1:n_inds){
  points(xs[i,seq(t0,tf,step)],ys[i,seq(t0,tf,step)],col=viridis(length(seq(t0,tf,step))),pch=19,cex=1)
  lines(xs[i,seq(t0,tf)],ys[i,seq(t0,tf)],col = ind_cols[i],lwd=1)
}

#Number of inds tracked over time
n_tracked <- colSums(!is.na(xs))
plot(timestamps, n_tracked)

dir.create(paste0('~/Desktop/lion_viz/day',day))
cocomo::generate_movement_and_calls_visualization(xs=xs, ys=ys, timestamps=timestamps, start_time = t0, end_time = tf, 
                                                  time_step = step, 
                                                  output_dir = paste0('~/Desktop/lion_viz/day',day),
                                                  tail_time=step*30,
                                                  calls=data.frame(), 
                                                  show_legend_calls = F, 
                                                  show_legend_inds = T,
                                                  ind_names = ids$code,
                                                  scalebar_size = 1000,
                                                  show_time = T)
