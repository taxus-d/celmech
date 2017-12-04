#!/usr/bin/gnuplot

if (exists("type")) {
        print type
} else {
        type = 'rk'
}

filen = 'data/'.type.'.dat'
outfile = 'img/'.type.'png'
set output outfile
set term pngcairo size 1000,1000

set grid
set xrange[-1.5:1.5]
set yrange[-1.5:1.5]

plot filen usi 2:3 w l
system('feh '.outfile.'&')
