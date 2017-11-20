#!/usr/bin/gnuplot

tim = 0
end_time = 1
system('mkdir -p pert')
# set samples 10000
set term png size 200, 200
#set term gif animate optimize
#set output 'anim3Dbody.gif'
#load 'spir.plt'
#set output 'animPertubations.gif'
load 'pertubations.plt'
#system('gwenview animPertubations.gif')
