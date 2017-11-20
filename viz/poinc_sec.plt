#!/usr/bin/gnuplot

if (!exists("type")) type = 'rk'
filen = 'data/'.type.'.dat'
outfile = 'img/'.type.'png'
set output outfile
set term pngcairo size 1000,1000

set grid
set xrange[0:120]
set yrange[-1:1]

plot filen usi 14:3 w l
system('feh '.outfile)
