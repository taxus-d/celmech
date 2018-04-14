#!/usr/bin/gnuplot

filename = 'data/orbit-orig.dat'
nrows    = system(sprintf("wc -l %s | awk '{print $1}'",filename))
gifdir   = 'img/proc'
# system('mkdir -p pert')
# set samples 10000
#set term qt
set term pngcairo size 500, 500
set xrange[-1.5:1.5]
set yrange[-1.5:1.5]

system("mkdir -p ".gifdir)
do for [ii=1:nrows] {
    framename = sprintf("%s/proc%05d.png",gifdir,ii)
    set output framename
     #'data/rk.dat' every ::1::ii using 2:3  with lines ls 1 notitle, 
    plot filename every ::ii::ii using 2:3 with points pointtype 7 pointsize 1 lc 1 title "1", \
         filename every ::ii::ii using 4:5 with points pointtype 7 pointsize 1 lc 2 title "2", \
         filename every ::ii::ii using 6:7 with points pointtype 7 pointsize 1 lc 3 title "3"
}
# system ("cd ".gifdir."; ffmpeg -r 50 -i 'proc%05d.png' video.mp4 ")

# vim:sw=4
