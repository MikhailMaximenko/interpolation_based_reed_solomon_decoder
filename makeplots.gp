set encoding cp1251
unset title
set xlabel "длина кода"
set ylabel "количество операций"
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Arial, 12"
set output 'comparison_with_berlecamp.png'
plot 'interpolation_results_0_5rate_max_errs_long.txt' using ($1):($4 + $5) with linespoints title 'улучшенный алгоритм Кадира-Линя-Рознеса' noenhanced, \
	 'berlecamp_results_0_5rate_max_errs_long.txt' using ($1):($4 + $5) with linespoints title 'алгоритм Берлекжмпа-Месси' noenhanced, \
