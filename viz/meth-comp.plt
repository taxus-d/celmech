#!/usr/bin/gnuplot

filen_rk = 'data/rk.dat'
filen_ae = 'data/ae.dat'
filen_ai = 'data/ai.dat'
outfile = 'img/comp.png'
set output outfile
set term pngcairo size 1000,1000

set xrange[-1.5:1.5]
set yrange[-1.5:1.5]
set size ratio -1

plot filen_rk usi 2:3 w l title "rk", filen_ae usi 2:3 w l title "ae",  filen_ai usi 2:3 w l title "ai"

system('feh '.outfile.'&')
