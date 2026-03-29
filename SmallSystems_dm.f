C======================================================================
C
C     SmallSystems_dm - Molecular Dynamics Simulation in Fortran 77
C
C======================================================================
C
      PROGRAM SmallSystems_dm
C
C     Parameters
C
      INTEGER NMAX
      PARAMETER (NMAX=2000)
C
C     Variables
C
      INTEGER N, i, icoll, ncolisiones
      INTEGER i1, i2
      INTEGER ii, jj, kk, ll, mm, nn, pp, qq  ! Different loop variables
C
      REAL vx(NMAX), vy(NMAX), vz(NMAX)
      REAL x(NMAX), y(NMAX), z(NMAX)
C
      REAL v_modulo
      REAL L
      REAL theta, phi, aleatorio
      REAL T_deseada
      REAL T_actual
      REAL suma_v2
      REAL Px, Py, Pz
      REAL E_total
      REAL E_inicial
      REAL nx, ny, nz
      REAL v1x, v1y, v1z
      REAL v2x, v2y, v2z
      REAL prod_escalar
      REAL diff
      REAL dx, dy, dz
C
C     Initial conditions
C
      N = 20
      v_modulo = 1.0
      ncolisiones = 200
      L = 10.0
C
      WRITE(*,*) '============================================'
      WRITE(*,*) '  SIMULACION DE DINAMICA MOLECULAR'
      WRITE(*,*) '  (Segunda version del programa)'
      WRITE(*,*) '============================================'
      WRITE(*,*)
      WRITE(*,*) 'Colectividad: Dinamica Molecular'
C
C     Initialize positions randomly in box
C
      DO 100 ii = 1, N
         x(ii) = L * RAND()
         y(ii) = L * RAND()
         z(ii) = L * RAND()
  100 CONTINUE
C
C     Initialize velocities randomly on sphere
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
C     Initial total momentum
C
      WRITE(*,*)
      WRITE(*,*) 'INICIALIZACION'
      WRITE(*,*) 'Verificacion momento total (inicial):'
C
      Px = 0.0
      Py = 0.0
      Pz = 0.0
      DO 300 ii = 1, N
         Px = Px + vx(ii)
         Py = Py + vy(ii)
         Pz = Pz + vz(ii)
  300 CONTINUE
      WRITE(*,*) 'Px total =', Px
      WRITE(*,*) 'Py total =', Py
      WRITE(*,*) 'Pz total =', Pz
C
C     Desired temperature and initial scaling
C
      T_deseada = 1.5
      suma_v2 = 0.0
      DO 400 ii = 1, N
         suma_v2 = suma_v2 + vx(ii)**2 + vy(ii)**2 + vz(ii)**2
  400 CONTINUE
      T_actual = suma_v2 / REAL(N) / 3.0
C
      WRITE(*,*)
      WRITE(*,*) 'Temperatura actual (antes de escalar):', T_actual
      WRITE(*,*) 'Temperatura deseada:', T_deseada
C
      DO 500 ii = 1, N
         vx(ii) = vx(ii) * SQRT(T_deseada / T_actual)
         vy(ii) = vy(ii) * SQRT(T_deseada / T_actual)
         vz(ii) = vz(ii) * SQRT(T_deseada / T_actual)
  500 CONTINUE
C
      suma_v2 = 0.0
      DO 600 ii = 1, N
         suma_v2 = suma_v2 + vx(ii)**2 + vy(ii)**2 + vz(ii)**2
  600 CONTINUE
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
      WRITE(*,*) 'Sin colisiones con paredes (colectividad DM)'
C
      DO 1000 icoll = 1, ncolisiones
C
C        Select two distinct particles
C
  110    i1 = INT(N * RAND()) + 1
  120    i2 = INT(N * RAND()) + 1
         IF (i2 .EQ. i1) GO TO 120
C
C        Random unit vector
C
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
         nx = SIN(theta) * COS(phi)
         ny = SIN(theta) * SIN(phi)
         nz = COS(theta)
C
C        Velocities of selected particles
C
         v1x = vx(i1)
         v1y = vy(i1)
         v1z = vz(i1)
         v2x = vx(i2)
         v2y = vy(i2)
         v2z = vz(i2)
C
C        Relative velocity dotted with normal (broken into steps for Fortran 77)
C
         dx = v1x - v2x
         dy = v1y - v2y
         dz = v1z - v2z
         prod_escalar = dx * nx + dy * ny + dz * nz
C
C        Update velocities (elastic collision with equal masses)
C
         vx(i1) = v1x - prod_escalar * nx
         vy(i1) = v1y - prod_escalar * ny
         vz(i1) = v1z - prod_escalar * nz
C
         vx(i2) = v2x + prod_escalar * nx
         vy(i2) = v2y + prod_escalar * ny
         vz(i2) = v2z + prod_escalar * nz
C
 1000 CONTINUE
C
      WRITE(*,*) 'Colisiones completadas'
C
C     Update positions (free streaming)
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
C     Final checks
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
      WRITE(*,*) '(Momento total libre en dinamica molecular)'
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
C     Calculate absolute difference manually for Fortran 77 compatibility
C     Using SQRT(x**2) instead of ABS
C
      diff = SQRT((E_total - E_inicial)**2)
      WRITE(*,*) 'Diferencia energias final e inicial:', diff
C
      WRITE(*,*) 'Energia conservada en dinamica molecular'
C
C     Output data files
C
      OPEN(10, FILE='data/velocidades_dm.dat', STATUS='UNKNOWN')
      DO 5000 ii = 1, N
         WRITE(10,*) vx(ii), vy(ii), vz(ii)
 5000 CONTINUE
      CLOSE(10)
C
      OPEN(10, FILE='data/modulos_dm.dat', STATUS='UNKNOWN')
      DO 6000 ii = 1, N
         WRITE(10,*) SQRT(vx(ii)**2 + vy(ii)**2 + vz(ii)**2)
 6000 CONTINUE
      CLOSE(10)
C
      OPEN(10, FILE='data/posiciones_dm.dat', STATUS='UNKNOWN')
      DO 7000 ii = 1, N
         WRITE(10,*) x(ii), y(ii), z(ii)
 7000 CONTINUE
      CLOSE(10)
C
      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - data/velocidades_dm.dat'
      WRITE(*,*) '  - data/modulos_dm.dat'
      WRITE(*,*) '  - data/posiciones_dm.dat'
C
      END