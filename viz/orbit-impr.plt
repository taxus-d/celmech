#!/usr/bin/gnuplot

if (exists("type")) {
        print type
} else {
        type = 'rk'
}
filen = 'data/'.type.'.dat'
filen1 = 'data/'.type.'1.dat'
outfile = 'img/'.type.'.png'
outfile1 = 'img/'.type.'1.png'
set term pngcairo size 1000,1000
set grid

set xrange[-1.5:1.5]
set yrange[-1.5:1.5]
set size ratio -1

set output outfile
plot filen usi 2:3 w l title "1", '' usi 4:5 w l title "2", '' usi 6:7 w l  title "3"

set output outfile1
plot filen1 usi 2:3 w l title "1", '' usi 4:5 w l title "2", '' usi 6:7 w l  title "3"

system('feh '.outfile.'&')
system('feh '.outfile1.'&')
