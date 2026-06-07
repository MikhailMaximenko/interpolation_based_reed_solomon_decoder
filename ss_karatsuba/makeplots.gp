set encoding utf8
unset title 
set xlabel 'Length'
set ylabel 'Operations'
unset grid
set logscale y
set logscale x
set xtics
set ytics
set terminal pngcairo size 1024,768 font "Arial, 12"
set output 'ss_karatsuba_mult.png'
plot 'schonhage_strassen_karatsuba9.txt' using ($1):($2 + $3) with linespoints title '3' noenhanced, \
    'schonhage_strassen_karatsuba27.txt' using ($1):($2 + $3) with linespoints title '27' noenhanced, \
    'schonhage_strassen_karatsuba81.txt' using ($1):($2 + $3) with linespoints title '81' noenhanced, \
    'schonhage_strassen_karatsuba243.txt' using ($1):($2 + $3) with linespoints title '243' noenhanced, \
    'schonhage_strassen_karatsuba729.txt' using ($1):($2 + $3) with linespoints title '729' noenhanced, \
    'schonhage_strassen_karatsuba2187.txt' using ($1):($2 + $3) with linespoints title '2187' noenhanced, \
    'schonhage_strassen_karatsuba6561.txt' using ($1):($2 + $3) with linespoints title '6561' noenhanced, \
