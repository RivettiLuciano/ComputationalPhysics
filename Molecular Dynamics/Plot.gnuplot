#! gnuplot
# cd 'D:\Documents\2017 2do semestre\Fiscom\practico 2\dinamica molecular\Git'
set terminal pngcairo size 550,430 enhanced font 'Verdana,10'
set output 'LoschmidParadox.png'

set yrange [-1.5:1.5]
set xrange [0.2:6.5]
set ylabel 'Entropy'
set xlabel 'Time'



plot 'LoschmidtParadox.txt' u 1:2 title ''

