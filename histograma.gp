set terminal png size 800,600
set output 'data/histograma.png'
set xlabel '|v|'
set ylabel 'P(|v|)'
set title 'Distribucion de rapideces - T=1.5'
set grid
set boxwidth 0.8 relative
plot 'data/histograma.dat' with boxes fill solid 0.7 title 'Simulacion', \
     'data/maxwell_teorica.dat' with lines linecolor rgb 'red' linewidth 2 title 'Maxwell-Boltzmann'
