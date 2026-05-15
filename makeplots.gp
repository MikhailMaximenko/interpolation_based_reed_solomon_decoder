unset title
set xlabel 'field size'
set ylabel 'operations'
unset grid
unset logscale y
unset logscale x
set xtics
set ytics
set terminal png size 1024,768 enhanced
set output 'results.png'
plot for [i=1:|ARGV|] ARGV[i] using ($1):($4 + $5) with linespoints title ARGV[i] noenhanced