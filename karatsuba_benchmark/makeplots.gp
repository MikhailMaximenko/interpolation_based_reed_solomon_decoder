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
set output 'karatsuba.png'
plot 'karatsuba2.txt' using ($1):($2 + $3) with linespoints title 'k=2' noenhanced, \
    'karatsuba4.txt' using ($1):($2 + $3) with linespoints title 'k=4' noenhanced, \
    'karatsuba8.txt' using ($1):($2 + $3) with linespoints title 'k=8' noenhanced, \
    'karatsuba16.txt' using ($1):($2 + $3) with linespoints title 'k=16' noenhanced, \
    'karatsuba32.txt' using ($1):($2 + $3) with linespoints title 'k=32' noenhanced, \
    'karatsuba64.txt' using ($1):($2 + $3) with linespoints title 'k=64' noenhanced, \
    'karatsuba128.txt' using ($1):($2 + $3) with linespoints title 'k=128' noenhanced, \
