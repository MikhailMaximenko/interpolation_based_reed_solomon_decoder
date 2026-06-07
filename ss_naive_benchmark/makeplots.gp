set encoding cp1251
unset title 
set xlabel 'длина'
set ylabel 'количество операций'
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Arial, 12"
set output 'ss_naive_mult.png'
plot 'schonhage_strassen3.txt' using ($1):($2 + $3) with linespoints title '3' noenhanced, \
	'schonhage_strassen9.txt' using ($1):($2 + $3) with linespoints title '9' noenhanced, \
	'schonhage_strassen27.txt' using ($1):($2 + $3) with linespoints title '27' noenhanced, \
	'schonhage_strassen81.txt' using ($1):($2 + $3) with linespoints title '81' noenhanced, \
	'schonhage_strassen243.txt' using ($1):($2 + $3) with linespoints title '243' noenhanced
