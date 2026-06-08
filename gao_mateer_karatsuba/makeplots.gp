set encoding cp1251
unset title 
set xlabel 'Длина'
set ylabel 'Количество операций'
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Arial, 12"
set output 'cmp_gao_karatsuba.png'
plot 'gao_mateer1.txt' using ($1):($2 + $3) with linespoints title 'Алгоритм Гао-Матира' noenhanced, \
    'karatsuba16.txt' using ($1):($2 + $3) with linespoints title 'Алгоритм Карацубы, 16' noenhanced, \
    'karatsuba32.txt' using ($1):($2 + $3) with linespoints title 'Алгоритм Карацубы, 32' noenhanced, \
