set encoding cp1251
unset title
set xlabel "Длина"
set ylabel "Количество операций"
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Times New Roman, 12"
set output 'total_result.png'
plot 'int_half_max_unimproved1.txt' using ($1):($4 + $5) with linespoints title 'Начальный алгоритм' noenhanced, \
	 'interpolation_decoder_final.txt' using ($1):($4 + $5) with linespoints title 'Итоговый алгоритм' noenhanced, \
