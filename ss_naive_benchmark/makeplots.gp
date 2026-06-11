set encoding cp1251
unset title 
set xlabel 'Длина'
set ylabel 'Количество операций'
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Times New Roman, 12"
set output 'ss_naive_mult.png'
plot 'schonhage_strassen3.txt' using ($1):($2 + $3) with linespoints title 'k=3' noenhanced, \
	'schonhage_strassen9.txt' using ($1):($2 + $3) with linespoints title 'k=9' noenhanced, \
	'schonhage_strassen27.txt' using ($1):($2 + $3) with linespoints title 'k=27' noenhanced, \
	'schonhage_strassen81.txt' using ($1):($2 + $3) with linespoints title 'k=81' noenhanced, \
	'schonhage_strassen243.txt' using ($1):($2 + $3) with linespoints title 'k=243' noenhanced
