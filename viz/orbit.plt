#!/usr/bin/gnuplot

if (exists("type")) {
} else {
        type = 'orbit-orig'
}
filen = 'data/'.type.'.dat'
outfile = 'img/'.type.'.png'
set output outfile
set term pngcairo size 1000,1000

set xrange[-1.5:1.5]
set yrange[-1.5:1.5]
set size ratio -1

plot filen using 2:3 with lines title "1", filen usi 4:5 w l title "2", filen usi 6:7 w l  title "3"

system('feh '.outfile.'&')
