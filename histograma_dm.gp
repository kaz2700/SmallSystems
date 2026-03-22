set terminal png size 800,600
set output 'data/histograma_dm.png'
set xlabel '|v|'
set ylabel 'P(|v|)'
set title 'DM: Distribucion de rapideces - T=1.5'
set grid
set boxwidth 0.8 relative
plot 'data/histograma_dm.dat' with boxes fill solid 0.7 title 'Simulacion DM', \
     'data/maxwell_teorica.dat' with lines linecolor rgb 'red' linewidth 2 title 'Maxwell-Boltzmann'
