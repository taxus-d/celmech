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

plot filen usi 1:2 w l, filen usi 1:3 w l
system('feh '.outfile.'&')
