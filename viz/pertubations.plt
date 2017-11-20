#set parametric
#set samples 10000
#set term png size 1280, 1080
outfile = sprintf('pert/pert%03.0f.png',1000*tim)
set output outfile
set xrange [-1:1]
set yrange [-1:1]
set size square

set trange[0:2*pi]

petals = 131
wave = tim
fx(a) = sin(petals*a)*cos(a)
fy(a) = sin((petals + wave)*a)*sin(a)

plot fx(t), fy(t) lt rgb "#1E90FF" lw 1 title "pertubations"

tim = tim + 0.001
if (tim < end_time) reread;
