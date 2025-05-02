#Make some plots of trajectories

load('~/Dropbox/lions_at_night/data_namibia/processed/lion_omukutu_xy_level1.RData')

#PARAMS
step <- 60
output_dir <- '~/Desktop/lion_viz/'

#dirs
output_dir_anim <- paste0(output_dir, 'anim/')
output_dir_static <- paste0(output_dir, 'static/')

n_days <- length(unique(lubridate::date(timestamps)))
for(day in 1:n_days){
  print(paste('day',day))
  t0 <- 24*60*60*day - 12*60*60 + 1
  tf <- t0 + 24*60*60
  
  #--- STATIC PLOT ---
  print('generating plot')
  n_inds <- nrow(xs)
  ind_cols <- c('red','orange','darkgreen','blue','purple','magenta')
  
  xmin <- min(xs[,t0:tf], na.rm=T)
  xmax <- max(xs[,t0:tf], na.rm=T)
  ymin <- min(ys[,t0:tf], na.rm=T)
  ymax <- max(ys[,t0:tf], na.rm=T)
  
  png(filename = paste0(output_dir_static, 'omukutu_day',day,'.png'), width = 800, height = 800, units = 'px')
  
  plot(NULL, xlim=c(xmin, xmax), ylim = c(ymin, ymax),asp=1,xlab='Easting (m)',ylab = 'Northing (m)')
  verts <- seq(round(xmin, digits = -3)-5000, round(xmax, digits = -3)+5000, 1000)
  horizs <- seq(round(ymin, digits = -3)-5000, round(ymax, digits = -3)+5000, 1000)
  abline(h=horizs, col = '#00000011',lwd=2)
  abline(v=verts, col = '#00000011',lwd=2)
  for(i in 1:n_inds){
    points(xs[i,seq(t0,tf,step)],ys[i,seq(t0,tf,step)],col=viridis(length(seq(t0,tf,step))),pch=19,cex=1)
    lines(xs[i,seq(t0,tf)],ys[i,seq(t0,tf)],col = ind_cols[i],lwd=1)
  }
  dev.off()
  
  #---ANIMATION---

  cocomo::generate_movement_and_calls_visualization(xs=xs, ys=ys, timestamps=timestamps, start_time = t0, end_time = tf, 
                                                    time_step = step, 
                                                    output_dir = '~/Desktop/lion_viz/anim/',
                                                    tail_time=step*30,
                                                    calls=data.frame(), 
                                                    show_legend_calls = F, 
                                                    show_legend_inds = T,
                                                    ind_names = ids$code,
                                                    scalebar_size = 1000,
                                                    show_time = T,
                                                    ind_point_size = 1.5)
  dir <- paste0(output_dir_anim,'seq_',t0,'-',tf,'/')
  system(command = paste('cd', dir))
  system(command = paste0('ffmpeg -framerate 40 -i %d.png -r 16 -c:v libx264 -pix_fmt yuvj420p omukutu_viz_night_',day,'.mp4'))
  system(command = paste0('mv omukutu_viz_night_',day,'.mp4 ../omukutu_viz_night_', day, '.mp4'))
}
