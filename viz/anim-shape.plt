#!/usr/bin/gnuplot

filename = 'data/orbit-orig.dat'
nrows    = system(sprintf("wc -l %s | awk '{print $1}'",filename))
gifdir   = 'img/proc-shape'
# system('mkdir -p pert')
# set samples 10000
#set term qt
set term pngcairo size 1000, 1000
set xrange[-1:1]
set yrange[-1:1]
set zrange[-1:1]
set view equal xyz
set view 57,280
set size ratio -1
set grid
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'

system("mkdir -p ".gifdir)

do for [ii=1:nrows] {
    framename = sprintf("%s/proc%05d.png",gifdir,ii)
    set output framename
    set multiplot layout 2,2
    # 3d 
    splot filename every ::ii::ii using 2:3:4 with points pointtype 7 pointsize 1 lc 1 title "trajectory",\
          filename every ::1::ii using 2:3:4 with lines dashtype 2 notitle
    #
    plot filename every ::ii::ii using 2:3 with points pointtype 7 pointsize 1 lc 1 title "x-y-trajectory",\
         filename every ::1::ii using 2:3 with lines dashtype 2 notitle
    #
    plot filename every ::ii::ii using 2:4 with points pointtype 7 pointsize 1 lc 1 title "x-z-trajectory",\
          filename every ::1::ii using 2:4 with lines dashtype 2 notitle
    #
    plot filename every ::ii::ii using 3:4 with points pointtype 7 pointsize 1 lc 1 title "y-z-trajectory",\
          filename every ::1::ii using 3:4 with lines dashtype 2 notitle
    #
    unset multiplot
}
# system ("cd ".gifdir."; ffmpeg -r 50 -i 'proc%05d.png' video.mp4 ")

# vim:sw=4
