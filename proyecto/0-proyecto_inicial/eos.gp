# eos: P vs rho
set encoding utf8

set xlabel "{/Symbol r}"
set ylabel "P"
set key l t

p 'eos.dat' u 1:4 w lp dt 2 pt 7 noti

set terminal pngcairo color
set output "eos.png"
rep
set terminal qt
