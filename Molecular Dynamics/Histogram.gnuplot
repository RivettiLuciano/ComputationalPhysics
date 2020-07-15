#! gnuplot
set terminal pngcairo size 550,430 enhanced font 'Verdana,10'
set output 'Histogram.png'

set style data histogram
set style histogram cluster gap 0
set style fill solid border -1
set boxwidth 0.9
set title 'Histogram'
set xtics border in scale 0,0 nomirror rotate by -45  autojustify
plot 'histogram.txt' using 2:xtic(int($1*10)%5 == 0 ? stringcolumn(1) : '') title ''