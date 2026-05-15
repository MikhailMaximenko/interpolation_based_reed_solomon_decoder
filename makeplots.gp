unset title
fname = "berlecamp_decoder_results.txt"
set xlabel '-'
set ylabel '-'
unset grid
set xtics
set ytics
set terminal png size 1024,768 enhanced
set output 'results.png'
plot fname using ($1):($4 + $5) with linespoints 