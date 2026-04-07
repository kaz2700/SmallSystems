C======================================================================
C
C     smallSystems - Molecular Dynamics Simulation in Fortran 77
C
C======================================================================
C
      PROGRAM smallSystems
C
C     Parameters
C
      INTEGER NMAX
      PARAMETER (NMAX=2000)
C
C     Variables
C
      INTEGER N, i, icoll, ncolisiones, icuantas
      INTEGER i1, i2, iwall, pared
      INTEGER ii, jj, kk
C
C     vx, vy, vz: componentes de velocidad de cada particula
C     x, y, z:     posiciones de cada particula
      REAL vx(NMAX), vy(NMAX), vz(NMAX)
      REAL x(NMAX), y(NMAX), z(NMAX)
C
C     Parametros del sistema
      REAL v_modulo         ! Velocidad inicial de cada particula
      REAL vcm_x, vcm_y, vcm_z  ! Velocidad del centro de masas
      REAL L               ! Lado del cubo (condiciones periodicas)
C
C     Variables para numeros aleatorios
      REAL theta, phi, aleatorio
C
C     Variables termodinamicas
      REAL T_deseada       ! Temperatura objetivo
      REAL T_actual        ! Temperatura actual (calculada)
      REAL suma_v2         ! Sumatoria de v^2
      REAL Px, Py, Pz      ! Componentes del momento total
      REAL E_total          ! Energia total
      REAL E_inicial        ! Energia inicial (para verificar conservacion)
C
C     Variables para colisiones elasticas
      REAL nx, ny, nz      ! Componentes del eje de colision (vector unitario)
      REAL v1x, v1y, v1z   ! Velocidad de particula 1 antes de colision
      REAL v2x, v2y, v2z   ! Velocidad de particula 2 antes de colision
      REAL prod_escalar     ! Producto escalar (v1-v2)·n
C
C     Variables para colisiones con pared
      REAL probabilidad     ! Probabilidad de colision con pared
C
C     Variables auxiliares
      REAL dx, dy, dz
      REAL diff
      REAL v_max
      REAL counts(100)
      REAL bin_width, bin_center
      INTEGER nbins, ibin
C
C     Condiciones iniciales
C
      N = 20                 ! Numero de particulas
      v_modulo = 1.0         ! Velocidad inicial (magnitud)
      ncolisiones = 200      ! Numero de colisiones a simular (10 por particula)
      icuantas = N           ! Cada N colisiones particula-particula, colisiones con pared
      probabilidad = 0.5    ! Probabilidad de colision con pared por particula
      L = 10.0               ! Tamano de la caja
C
C     Inicializar posiciones aleatoriamente en la caja
C
      DO 100 ii = 1, N
         x(ii) = L * RAND()
         y(ii) = L * RAND()
         z(ii) = L * RAND()
  100 CONTINUE
C
C     Para una distribucion uniforme sobre la esfera unitaria, usamos:
C     - theta = arccos(2*u - 1) con u~U(0,1)  -> distribucion uniforme en [0, pi]
C     - phi = 2*pi*u con u~U(0,1)             -> distribucion uniforme en [0, 2pi]
C
      DO 200 ii = 1, N
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
C
         vx(ii) = v_modulo * SIN(theta) * COS(phi)
         vy(ii) = v_modulo * SIN(theta) * SIN(phi)
         vz(ii) = v_modulo * COS(theta)
  200 CONTINUE
C
C     En el ensemble microcanonico, el momento total debe ser cero.
C     v_i' = v_i - V_cm
C     donde V_cm = (1/N) * sum(v_i)
C
      vcm_x = 0.0
      vcm_y = 0.0
      vcm_z = 0.0
      DO 300 ii = 1, N
         vcm_x = vcm_x + vx(ii)
         vcm_y = vcm_y + vy(ii)
         vcm_z = vcm_z + vz(ii)
  300 CONTINUE
      vcm_x = vcm_x / REAL(N)
      vcm_y = vcm_y / REAL(N)
      vcm_z = vcm_z / REAL(N)
C
      DO 400 ii = 1, N
         vx(ii) = vx(ii) - vcm_x
         vy(ii) = vy(ii) - vcm_y
         vz(ii) = vz(ii) - vcm_z
  400 CONTINUE
C
C     Escalado de temperatura (T = <v^2> / 3)
C
      T_deseada = 1.5
      suma_v2 = 0.0
      DO 500 ii = 1, N
         suma_v2 = suma_v2 + vx(ii)**2 + vy(ii)**2 + vz(ii)**2
  500 CONTINUE
      T_actual = suma_v2 / REAL(N) / 3.0
C
      WRITE(*,*)
      WRITE(*,*) 'Temperatura actual (antes de escalar):', T_actual
      WRITE(*,*) 'Temperatura deseada:', T_deseada
C
      DO 600 ii = 1, N
         vx(ii) = vx(ii) * SQRT(T_deseada / T_actual)
         vy(ii) = vy(ii) * SQRT(T_deseada / T_actual)
         vz(ii) = vz(ii) * SQRT(T_deseada / T_actual)
  600 CONTINUE
C
C     Verificar energia inicial
C
      suma_v2 = 0.0
      DO 700 ii = 1, N
         suma_v2 = suma_v2 + vx(ii)**2 + vy(ii)**2 + vz(ii)**2
  700 CONTINUE
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_inicial = 0.5 * suma_v2
C
      WRITE(*,*)
      WRITE(*,*) 'Temperatura despues de escalar:', T_actual
      WRITE(*,*) 'Energia cinetica inicial:', E_inicial
C
C     Bucle de colisiones
C
      WRITE(*,*)
      WRITE(*,*) 'SIMULANDO COLISIONES'
      WRITE(*,*) 'Numero de colisiones particula-particula:', ncolisiones
      WRITE(*,*) 'Colisiones con pared cada:', icuantas
C
      DO 1000 icoll = 1, ncolisiones
C
C     Colision particula-particula
C     Seleccionar dos particulas distintas
C
  110    i1 = INT(N * RAND()) + 1
  120    i2 = INT(N * RAND()) + 1
         IF (i2 .EQ. i1) GO TO 120
C
C     Generar eje de colision aleatorio
C
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
         nx = SIN(theta) * COS(phi)
         ny = SIN(theta) * SIN(phi)
         nz = COS(theta)
C
C     Guardar velocidades antes de la colision
C
         v1x = vx(i1)
         v1y = vy(i1)
         v1z = vz(i1)
         v2x = vx(i2)
         v2y = vy(i2)
         v2z = vz(i2)
C
C     Aplicar formula de colision elastica (dividida para Fortran 77)
C
         dx = v1x - v2x
         dy = v1y - v2y
         dz = v1z - v2z
         prod_escalar = dx * nx + dy * ny + dz * nz
C
         vx(i1) = v1x - prod_escalar * nx
         vy(i1) = v1y - prod_escalar * ny
         vz(i1) = v1z - prod_escalar * nz
C
         vx(i2) = v2x + prod_escalar * nx
         vy(i2) = v2y + prod_escalar * ny
         vz(i2) = v2z + prod_escalar * nz
C
C     Colisiones con pared: cada icuantas colisiones, seleccionar N particulas
C     y cada una tiene probabilidad de chocar con una pared
C
         IF (MOD(icoll, icuantas) .EQ. 0) THEN
            DO 900 iwall = 1, N
               aleatorio = RAND()
               IF (aleatorio .LT. probabilidad) THEN
                  i1 = INT(N * RAND()) + 1
                  pared = INT(6.0 * RAND()) + 1
                  IF (pared .EQ. 1) vx(i1) = -vx(i1)
                  IF (pared .EQ. 2) vx(i1) = -vx(i1)
                  IF (pared .EQ. 3) vy(i1) = -vy(i1)
                  IF (pared .EQ. 4) vy(i1) = -vy(i1)
                  IF (pared .EQ. 5) vz(i1) = -vz(i1)
                  IF (pared .EQ. 6) vz(i1) = -vz(i1)
               END IF
  900       CONTINUE
         END IF
C
 1000 CONTINUE
C
      WRITE(*,*) 'Colisiones completadas'
C
C     Actualizar posiciones
C
      WRITE(*,*)
      WRITE(*,*) 'ACTUALIZANDO POSICIONES'
C
      DO 2000 ii = 1, N
         x(ii) = x(ii) + vx(ii)
         y(ii) = y(ii) + vy(ii)
         z(ii) = z(ii) + vz(ii)
 2000 CONTINUE
C
      WRITE(*,*) 'Posiciones actualizadas'
C
C     Verificar resultados
C
      WRITE(*,*)
      WRITE(*,*) 'VERIFICACION RESULTADOS'
C
      Px = 0.0
      Py = 0.0
      Pz = 0.0
      DO 3000 ii = 1, N
         Px = Px + vx(ii)
         Py = Py + vy(ii)
         Pz = Pz + vz(ii)
 3000 CONTINUE
      WRITE(*,*) 'Momento total final:'
      WRITE(*,*) 'Px =', Px, '  Py =', Py, '  Pz =', Pz
C
C     Verificar conservacion de energia
C
      suma_v2 = 0.0
      DO 4000 ii = 1, N
         suma_v2 = suma_v2 + vx(ii)**2 + vy(ii)**2 + vz(ii)**2
 4000 CONTINUE
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_total = 0.5 * suma_v2
C
      WRITE(*,*)
      WRITE(*,*) 'Temperatura final:', T_actual
      WRITE(*,*) 'Energia cinetica final:', E_total
C
C     Valor absoluto manual (compatible con Fortran 77)
C
      diff = E_total - E_inicial
      IF (diff .LT. 0.0) THEN
         diff = -diff
      END IF
      WRITE(*,*) 'Diferencia energias final e inicial:', diff
      WRITE(*,*)
      WRITE(*,*) 'Generando histograma...'
C
C     Calcular velocidad maxima y ancho de bins
C
      v_modulo = SQRT(8.0 * T_deseada / 3.141592653589793)
      v_max = 0.0
      DO 7100 ii = 1, N
         v_modulo = SQRT(vx(ii)**2 + vy(ii)**2 + vz(ii)**2)
         IF (v_modulo .GT. v_max) v_max = v_modulo
 7100 CONTINUE
      v_max = v_max * 1.2
      nbins = 30
      bin_width = v_max / REAL(nbins)
C
C     Inicializar contadores de bins a cero
C
      DO 7200 ii = 1, nbins
         counts(ii) = 0.0
 7200 CONTINUE
C
C     Contar particulas en cada bin
C
      DO 7300 ii = 1, N
         v_modulo = SQRT(vx(ii)**2 + vy(ii)**2 + vz(ii)**2)
         ibin = INT(v_modulo / bin_width) + 1
         IF (ibin .GT. nbins) ibin = nbins
         IF (ibin .LT. 1) ibin = 1
         counts(ibin) = counts(ibin) + 1.0
 7300 CONTINUE
C
C     Normalizar histograma
C
      DO 7400 ii = 1, nbins
         counts(ii) = counts(ii) / (REAL(N) * bin_width)
 7400 CONTINUE
C
C     Escribir histograma a archivo
C
      OPEN(10, FILE='data/histograma.dat', STATUS='UNKNOWN')
      DO 7500 ii = 1, nbins
         bin_center = (REAL(ii) - 0.5) * bin_width
         WRITE(10,8010) bin_center, counts(ii)
 7500 CONTINUE
      CLOSE(10)
 8010 FORMAT(F7.4, 1X, F10.6)

      WRITE(*,*)
      WRITE(*,*) 'Generado: data/histograma.dat'

      END