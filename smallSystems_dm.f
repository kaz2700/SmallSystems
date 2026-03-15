      PROGRAM smallSystems
      IMPLICIT NONE

!======================================================================
!  PROGRAMA v2: DINAMICA MOLECULAR
!  Simulacion de N particulas con integracion temporal
!  Condiciones de contorno periodicas (PBC)
!  Autor: 
!  Fecha: 
!======================================================================

      INTEGER, PARAMETER :: NMAX = 1000
      INTEGER :: N, i, j, icoll, nsteps, ncolisiones
      INTEGER :: i1, i2

!---------------------------------------------------------------------
!  Variables cinematicas
!---------------------------------------------------------------------
      REAL :: vx(NMAX), vy(NMAX), vz(NMAX)
      REAL :: x(NMAX), y(NMAX), z(NMAX)
      REAL :: vx_new(NMAX), vy_new(NMAX), vz_new(NMAX)
      REAL :: x_new(NMAX), y_new(NMAX), z_new(NMAX)
      REAL :: v_modulo, vcm_x, vcm_y, vcm_z
      REAL :: L, dt

!---------------------------------------------------------------------
!  Variables para angulos y numeros aleatorios
!---------------------------------------------------------------------
      REAL :: theta, phi, aleatorio

!---------------------------------------------------------------------
!  Variables termodinamicas
!---------------------------------------------------------------------
      REAL :: T_deseada, T_actual, suma_v2
      REAL :: Px, Py, Pz, E_total, E_inicial, E_pot
      REAL :: E_inicial_tot, E_final_tot

!---------------------------------------------------------------------
!  Variables para colisiones (opcional para gas con interaccion)
!---------------------------------------------------------------------
      REAL :: nx, ny, nz
      REAL :: v1x, v1y, v1z, v2x, v2y, v2z
      REAL :: prod_escalar

!======================================================================
!  SECCION 3.1: INICIALIZACION
!======================================================================

!  Parametros del sistema
      N = 100
      L = 10.0
      dt = 0.01
      nsteps = 1000
      ncolisiones = 0

!---------------------------------------------------------------------
!  3.1.1: Inicializar posiciones aleatorias en el volumen
!---------------------------------------------------------------------
      DO i = 1, N
         x(i) = L * RAND()
         y(i) = L * RAND()
         z(i) = L * RAND()
      END DO

!---------------------------------------------------------------------
!  3.1.2: Inicializar velocidades con igual modulo y direccion aleatoria
!---------------------------------------------------------------------
      v_modulo = 1.0
      DO i = 1, N
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio

         vx(i) = v_modulo * SIN(theta) * COS(phi)
         vy(i) = v_modulo * SIN(theta) * SIN(phi)
         vz(i) = v_modulo * COS(theta)
      END DO

!---------------------------------------------------------------------
!  3.1.3: Imponer momento total cero
!---------------------------------------------------------------------
      vcm_x = SUM(vx(1:N)) / REAL(N)
      vcm_y = SUM(vy(1:N)) / REAL(N)
      vcm_z = SUM(vz(1:N)) / REAL(N)

      DO i = 1, N
         vx(i) = vx(i) - vcm_x
         vy(i) = vy(i) - vcm_y
         vz(i) = vz(i) - vcm_z
      END DO

!---------------------------------------------------------------------
!  3.1.4: Escalar velocidades para alcanzar temperatura deseada
!---------------------------------------------------------------------
      T_deseada = 1.5
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))

      DO i = 1, N
         vx(i) = vx(i) * SQRT(T_deseada / T_actual)
         vy(i) = vy(i) * SQRT(T_deseada / T_actual)
         vz(i) = vz(i) * SQRT(T_deseada / T_actual)
      END DO

!---------------------------------------------------------------------
!  Calcular energia inicial
!---------------------------------------------------------------------
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      E_inicial = 0.5 * suma_v2
      E_pot = 0.0
      E_inicial_tot = E_inicial + E_pot

      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))

      WRITE(*,*) '=== INICIALIZACION ==='
      WRITE(*,*) 'N =', N, '  L =', L, '  dt =', dt
      WRITE(*,*) 'nsteps =', nsteps
      WRITE(*,*) 'Temperatura inicial:', T_deseada
      WRITE(*,*) 'Energia cinetica inicial:', E_inicial
      WRITE(*,*) 'Momento total inicial: Px =', Px, ' Py =', Py, ' Pz =', Pz

!======================================================================
!  SECCION 3.2: DINAMICA MOLECULAR - INTEGRACION TEMPORAL
!======================================================================

!---------------------------------------------------------------------
!  Algoritmo de Velocity Verlet:
!  1. r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
!  2. v(t+dt/2) = v(t) + 0.5*a(t)*dt
!  3. calcular fuerzas a(t+dt)
!  4. v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
!
!  Para gas ideal sin interacciones: a = 0
!  Las particulas se mueven en linha recta con PBC
!---------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,*) '=== DINAMICA MOLECULAR ==='

      DO icoll = 1, nsteps
!        Actualizar posiciones: x(t+dt) = x(t) + v*dt
         DO i = 1, N
            x(i) = x(i) + vx(i) * dt
            y(i) = y(i) + vy(i) * dt
            z(i) = z(i) + vz(i) * dt
         END DO

!        Aplicar PBC
         DO i = 1, N
            IF (x(i) > L) x(i) = x(i) - L
            IF (x(i) < 0.0) x(i) = x(i) + L
            IF (y(i) > L) y(i) = y(i) - L
            IF (y(i) < 0.0) y(i) = y(i) + L
            IF (z(i) > L) z(i) = z(i) - L
            IF (z(i) < 0.0) z(i) = z(i) + L
         END DO
      END DO

      WRITE(*,*) 'Integracion completada:', nsteps, 'pasos'

!======================================================================
!  SECCION 3.3: VERIFICACION
!======================================================================

!---------------------------------------------------------------------
!  3.2.1: Verificar conservacion de momento total
!---------------------------------------------------------------------
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))

      WRITE(*,*)
      WRITE(*,*) '=== VERIFICACION ==='
      WRITE(*,*) 'Momento total final:'
      WRITE(*,*) 'Px =', Px, '  Py =', Py, '  Pz =', Pz

!---------------------------------------------------------------------
!  3.2.2: Verificar conservacion de energia
!---------------------------------------------------------------------
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      E_total = 0.5 * suma_v2
      E_pot = 0.0
      E_final_tot = E_total + E_pot

      WRITE(*,*)
      WRITE(*,*) 'Energia cinetica final:', E_total
      WRITE(*,*) 'Energia total final:', E_final_tot
      WRITE(*,*) 'Diferencia energia:', ABS(E_final_tot - E_inicial_tot)

!---------------------------------------------------------------------
!  3.2.3: Comparar distribucion con teoria microcanonica
!---------------------------------------------------------------------
      OPEN(10, FILE='velocidades_dm.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) vx(i), vy(i), vz(i)
      END DO
      CLOSE(10)

      OPEN(10, FILE='modulos_dm.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) SQRT(vx(i)**2 + vy(i)**2 + vz(i)**2)
      END DO
      CLOSE(10)

      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - velocidades_dm.dat'
      WRITE(*,*) '  - modulos_dm.dat'

      END PROGRAM
