set parametric
#set term png size 1280, 1080
#outfile = sprintf('animation/spiral%03.0f.png',100*tim)
#set output outfile
set samples 1000
set xrange [-1:1]
set yrange [-1:1]
set zrange [-1:1]
val = 0.1*tim*18*pi
val0 = 0.1*18*pi
sp_size = 0.05
set urange [0:val]
set isosamples 10

fx(a) = sin(a/2)*cos(a)
fy(a) = sin(a/4)*sin(a)
fz(a) = sin(a/8)

splot fx(u), fy(u), fz(u) lt rgb "#1E90FF" lw 3 title "trajectory", \
sp_size*cos(u)*cos(v)+fx(val), sp_size*sin(u)*cos(v)+fy(val), sp_size*sin(v)+fz(val) lt rgb "#1E90FF" notitle

tim = tim + 0.03
if (tim < end_time) reread;
