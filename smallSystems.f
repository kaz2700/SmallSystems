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
C     Auxiliary variables
      REAL dx, dy, dz
      REAL diff
C
C     Initial conditions
C
      N = 20                 ! Numero de particulas
      v_modulo = 1.0         ! Velocidad inicial (magnitud)
      ncolisiones = 200      ! Numero de colisiones a simular (10 por particula)
      icuantas = N           ! Cada N colisiones particula-particula, colisiones con pared
      probabilidad = 0.5    ! Probabilidad de colision con pared por particula
      L = 10.0               ! Tamano de la caja
C
C     Initialize positions randomly in box
C
      DO 100 ii = 1, N
         x(ii) = L * RAND()
         y(ii) = L * RAND()
         z(ii) = L * RAND()
  100 CONTINUE
C
C     For a uniform distribution on the unit sphere, we use:
C     - theta = arccos(2*u - 1) with u~U(0,1)  -> uniform distribution in [0, pi]
C     - phi = 2*pi*u with u~U(0,1)             -> uniform distribution in [0, 2pi]
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
C     In the microcanonical ensemble, the total momentum must be zero.
C     v_i' = v_i - V_cm
C     where V_cm = (1/N) * sum(v_i)
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
C     Temperature scaling (T = <v^2> / 3)
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
C     Verify initial energy
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
C     Collision loop
C
      WRITE(*,*)
      WRITE(*,*) 'SIMULANDO COLISIONES'
      WRITE(*,*) 'Numero de colisiones particula-particula:', ncolisiones
      WRITE(*,*) 'Colisiones con pared cada:', icuantas
C
      DO 1000 icoll = 1, ncolisiones
C
C     Particle-particle collision
C     Select two distinct particles
C
  110    i1 = INT(N * RAND()) + 1
  120    i2 = INT(N * RAND()) + 1
         IF (i2 .EQ. i1) GO TO 120
C
C     Generate random collision axis
C
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
         nx = SIN(theta) * COS(phi)
         ny = SIN(theta) * SIN(phi)
         nz = COS(theta)
C
C     Save velocities before collision
C
         v1x = vx(i1)
         v1y = vy(i1)
         v1z = vz(i1)
         v2x = vx(i2)
         v2y = vy(i2)
         v2z = vz(i2)
C
C     Apply elastic collision formula (split for Fortran 77)
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
C     Wall collisions: every icuantas collisions, select N particles
C     and each has probability of colliding with a wall
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
C     Update positions
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
C     Verify results
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
C     Verify energy conservation
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
C     Manual absolute value (Fortran 77 compatible)
C
      diff = E_total - E_inicial
      IF (diff .LT. 0.0) THEN
         diff = -diff
      END IF
      WRITE(*,*) 'Diferencia energias final e inicial:', diff
C
C     Output data files
C
      OPEN(10, FILE='data/velocidades_final.dat', STATUS='UNKNOWN')
      DO 5000 ii = 1, N
         WRITE(10,*) vx(ii), vy(ii), vz(ii)
 5000 CONTINUE
      CLOSE(10)
C
      OPEN(10, FILE='data/modulos_final.dat', STATUS='UNKNOWN')
      DO 6000 ii = 1, N
         WRITE(10,*) SQRT(vx(ii)**2 + vy(ii)**2 + vz(ii)**2)
 6000 CONTINUE
      CLOSE(10)
C
      OPEN(10, FILE='data/posiciones_final.dat', STATUS='UNKNOWN')
      DO 7000 ii = 1, N
         WRITE(10,*) x(ii), y(ii), z(ii)
 7000 CONTINUE
      CLOSE(10)
C
      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - data/velocidades_final.dat (vx, vy, vz por linea)'
      WRITE(*,*) '  - data/modulos_final.dat (|v| por linea)'
      WRITE(*,*) '  - data/posiciones_final.dat (x, y, z por linea)'
C
      END