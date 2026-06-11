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
set output 'interpolation_karatsuba.png'
plot 'interpolation_half_max256.txt' using ($1):($4 + $5) with linespoints title 'k=256' noenhanced, \
    'interpolation_half_max512.txt' using ($1):($4 + $5) with linespoints title 'k=512' noenhanced, \
    'interpolation_half_max1024.txt' using ($1):($4 + $5) with linespoints title 'k=1024' noenhanced, \
    'interpolation_half_max2048.txt' using ($1):($4 + $5) with linespoints title 'k=2048' noenhanced, \
    'interpolation_half_max4096.txt' using ($1):($4 + $5) with linespoints title 'k=4096' noenhanced, \
