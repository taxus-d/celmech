if (exists("type")) {
} else {
        type = 'orbit-orig'
}
filen = 'data/'.type.'.dat'
outfile = 'img/'.type.'-shape'.'.png'
print filen
set output outfile
set term pngcairo size 1000,1000


set parametric
set urange[0:2*pi]
set vrange[-pi/2:pi/2]
r = 1
fx(u,v) = r*cos(v)*cos(u)
fy(u,v) = r*cos(v)*sin(u)
fz(u,v) = r*sin(v)
set isosamples 30
set view equal xyz

set style line 9 linetype rgb "gray" linewidth 1
set style line 1 linetype rgb "skyblue" linewidth 2
set style line 2 linetype rgb "steelblue" linewidth 2

splot fx(u,v),fy(u,v),fz(u,v) with lines linestyle 9,\
      filen using 2:3:4 with lines linestyle 1, \
      'data/orbit-ideal.dat' using 2:3:4 with lines linestyle 2

system('feh '.outfile.'&')
