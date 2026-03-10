Física Estadística Avanzada
Trabajo sobre Teoría Cinética de sistemas pequeños
1
Trabajo sobre Teoría Cinética de
sistemas pequeños
Teoría:
1.- Colectividad microcanónica.
A partir de la función densidad en el equilibrio en la colectividad microcanónica
Integrando sobre todos los grados de libertad menos los correspondientes al momento de una
partícula:
hay que obtener la función densidad para el momento de una partícula:
Tras realizar los cambios de coordenadas adecuados también hay que obtener la distribución
de velocidades y de energía.
Finalmente, tomando el límite de N tendiendo a infinito hay que recuperar el valor de la
distribución de velocidades de Maxwell-Boltzmann.
Física Estadística Avanzada
Trabajo sobre Teoría Cinética de sistemas pequeños
2
2.- Colectividad de dinámica molecular.
En este caso la función densidad es
Donde
Admitiendo el resultado de Román et al
hay que obtener la función densidad para el momento de una partícula en el caso de que el
momento total sea igual a cero
Física Estadística Avanzada
Trabajo sobre Teoría Cinética de sistemas pequeños
3
3.- Método de la transformada de Laplace.
La idea de este método es la siguiente:
A la hora de realizar cálculos en estas funciones se utiliza la técnica de la transformada de
Laplace que consiste en tomar la transformada de Laplace con respecto de la energía y en ese
caso se pueden hacer las integrales sobre los momentos de las partículas, si a continuación se
hace la transformada inversa se obtiene: (en el caso de la colectividad de dinámica molecular
también hay que tomar transformadas de fourier respecto del momento total del sistema)
Física Estadística Avanzada
Trabajo sobre Teoría Cinética de sistemas pequeños
4
En este apartado hay que:
• Demostrar las ecuaciones anteriores aplicando la técnica.
• Repetir los resultados para las distribuciones de velocidad y energía en la colectividad
microcanónica y de dinámica molecular por el método de la transformada de Laplace.
Práctica:
1.- Primera versión del programa (Microcanónica)
1. Escribir un programa que almacene las tres componentes de la velocidad de cada
particula en 3 arrays (con un máximo de 32 elementos, por ejemplo): vx(i),vy(i),vz(i)
2. Asignar velocidades iniciales con igual módulo y dirección aleatoria
3. Imponer momento total igual a cero
4. Escalar las velocidades para fijar la temperatura del sistema
5. Elegir pares de partículas y hacerlas colisionar como partículas elásticas
6. Cada cierto número de colisiones entre partículas (N) elegir N partículas y hacerlas
colisionar con una cierta probabilidad con una pared de tal manera que una de las
componentes de la velocidad cambie de signo.
7. Al cabo de un número suficientemente grande de colisiones (10 por partícula, por
ejemplo) muestrear la distribución de velocidades. Repetir el muestreo
acumulativamente hasta obtener un histograma ‘suficientemente suave’
8. Comparar con los resultados teóricos.
2.- Segunda versión del programa (Dinámica Molecular)
1. Escribir un programa que almacene las tres componentes de la velocidad de cada
particula en 3 arrays (con un máximo de 32 elementos, por ejemplo): vx(i),vy(i),vz(i)
2. Asignar velocidades iniciales con igual módulo y dirección aleatoria
3. Imponer momento total igual a cero
4. Escalar las velocidades para fijar la temperatura del sistema
5. Elegir pares de partículas y hacerlas colisionar como partículas elásticas
6. En este caso no hay que colisionar las partículas con las paredes del recinto, se puede
adaptar el programa anterior para que nos permita elegir si hay conservación de
momento o no
Física Estadística Avanzada
Trabajo sobre Teoría Cinética de sistemas pequeños
5
7. Al cabo de un número suficientemente grande de colisiones (10 por partícula, por
ejemplo) muestrear la distribución de velocidades. Repetir el muestreo
acumulativamente hasta obtener un histograma ‘suficientemente suave’
8. Comparar con los resultados teóricos.
3.- Tercera versión del programa
1. Introducir las modificaciones para que se ajuste a la descripción de DSMC tal y como
viene en el paper de Montanero y Santos. Comparar con los resultados de los
programas anteriores.
2. Modificar el programa para que muestree un gran número de sistemas
simultáneamente.
3. Comparar con la teoría en las dos colectividades consideradas previamente.
