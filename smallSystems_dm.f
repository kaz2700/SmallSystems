      PROGRAM smallSystems
      IMPLICIT NONE

!======================================================================
!  PROGRAMA v2: DINAMICA MOLECULAR (MD)
!  =====================================
!
!  Este programa simula un gas ideal usando Dinamica Molecular.
!  A diferencia del metodo de Monte Carlo (v1), aqui evolucionamos
!  el sistema explicitamente en el tiempo.
!
!  Caracteristicas principales:
!  1. Integracion temporal usando algoritmo de Verlet
!  2. Sin fuerzas externas (gas ideal): a = 0
!  3. Las particulas se mueven en linea recta hasta chocar
!  4. Condiciones de contorno periodicas (PBC)
!
!  DIFERENCIA CON v1 (Microcanonica):
!  ----------------------------------
!  - v1: Colisiones aleatorias discretas, evolucion no temporal
!  - v2: Evolucion temporal explícita con dt
!
!  En un gas ideal sin interacciones, las particulas NUNCA cambian
!  su velocidad. Por lo tanto, la distribucion NO sera Maxwell-
!  Boltzmann (solo lo seria si hubiera colisiones que redistribuyan
!  la energia).
!
!  Para ejecutar:
!    gfortran -ffree-form smallSystems_dm.f -o smallSystems_dm
!    ./smallSystems_dm
!
!======================================================================

      INTEGER, PARAMETER :: NMAX = 1000
      INTEGER :: N, i, istep, nsteps
      INTEGER :: i1, i2

!---------------------------------------------------------------------
!  DECLARACION DE VARIABLES
!---------------------------------------------------------------------
!  vx, vy, vz: componentes de velocidad de cada particula
!  x, y, z:     posiciones de cada particula
      REAL :: vx(NMAX), vy(NMAX), vz(NMAX)
      REAL :: x(NMAX), y(NMAX), z(NMAX)

!  Parametros del sistema
      REAL :: v_modulo         ! Velocidad inicial
      REAL :: vcm_x, vcm_y, vcm_z  ! Velocidad del centro de masas
      REAL :: L               ! Tamano de la caja
      REAL :: dt              ! Paso de tiempo

!  Variables para numeros aleatorios
      REAL :: theta, phi, aleatorio

!  Variables termodinamicas
      REAL :: T_deseada       ! Temperatura objetivo
      REAL :: T_actual        ! Temperatura actual
      REAL :: suma_v2         ! Sumatoria de v^2
      REAL :: Px, Py, Pz      ! Momento total
      REAL :: E_cinetica      ! Energia cinetica
      REAL :: E_pot           ! Energia potencial (0 para gas ideal)
      REAL :: E_total         ! Energia total

!======================================================================
!  SECCION 3.1: INICIALIZACION
!======================================================================
!  Configuracion inicial del sistema:
!  1. Definir parametros (N, L, dt, numero de pasos)
!  2. Asignar posiciones aleatorias
!  3. Asignar velocidades aleatorias
!  4. Imponer momento total cero
!  5. Escalar a temperatura deseada
!======================================================================

!  --- Parametros del sistema ---
      N = 100                 ! Numero de particulas
      L = 10.0               ! Tamano de la caja
      dt = 0.01              ! Paso de tiempo
      nsteps = 1000          ! Numero de pasos de integracion

!  --- 3.1.1: Inicializar posiciones aleatorias ---
      DO i = 1, N
         x(i) = L * RAND()
         y(i) = L * RAND()
         z(i) = L * RAND()
      END DO

!  --- 3.1.2: Inicializar velocidades aleatorias ---
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

!  --- 3.1.3: Imponer momento total cero ---
      vcm_x = SUM(vx(1:N)) / REAL(N)
      vcm_y = SUM(vy(1:N)) / REAL(N)
      vcm_z = SUM(vz(1:N)) / REAL(N)

      DO i = 1, N
         vx(i) = vx(i) - vcm_x
         vy(i) = vy(i) - vcm_y
         vz(i) = vz(i) - vcm_z
      END DO

!  --- 3.1.4: Escalar velocidades a temperatura deseada ---
      T_deseada = 1.5
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))

      DO i = 1, N
         vx(i) = vx(i) * SQRT(T_deseada / T_actual)
         vy(i) = vy(i) * SQRT(T_deseada / T_actual)
         vz(i) = vz(i) * SQRT(T_deseada / T_actual)
      END DO

!  Calcular energia inicial
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      E_cinetica = 0.5 * suma_v2
      E_pot = 0.0
      E_total = E_cinetica + E_pot

      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))

      WRITE(*,*) '========================================'
      WRITE(*,*) '  PROGRAMA v2: DINAMICA MOLECULAR'
      WRITE(*,*) '========================================'
      WRITE(*,*)
      WRITE(*,*) '=== PARAMETROS ==='
      WRITE(*,*) 'Numero de particulas: N =', N
      WRITE(*,*) 'Tamano de caja:       L =', L
      WRITE(*,*) 'Paso de tiempo:       dt =', dt
      WRITE(*,*) 'Pasos de integracion: nsteps =', nsteps
      WRITE(*,*)
      WRITE(*,*) '=== ESTADO INICIAL ==='
      WRITE(*,*) 'Temperatura inicial:', T_deseada
      WRITE(*,*) 'Energia cinetica:', E_cinetica
      WRITE(*,*) 'Momento total: Px =', Px, ' Py =', Py, ' Pz =', Pz
      WRITE(*,*)
      WRITE(*,*) 'NOTA: En un gas ideal sin interacciones,'
      WRITE(*,*) 'las velocidades NO cambian con el tiempo.'
      WRITE(*,*) 'La distribucion no sera Maxwell-Boltzmann.'

!======================================================================
!  SECCION 3.2: INTEGRACION TEMPORAL
!======================================================================
!
!  ALGORITMO DE VELOCITY VERLET (para referencia):
!  -----------------------------------------------
!  1. r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
!  2. v(t+dt/2) = v(t) + 0.5*a(t)*dt
!  3. Calcular nuevas fuerzas a(t+dt)
!  4. v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
!
!  PARA GAS IDEAL (sin interacciones, a = 0):
!  -------------------------------------------
!  Las particulas se mueven en linea recta:
!    x(t+dt) = x(t) + v(t) * dt
!
!  Las velocidades permanecen CONSTANTES durante toda la simulacion.
!
!======================================================================

      WRITE(*,*)
      WRITE(*,*) '=== EVOLUCION TEMPORAL ==='
      WRITE(*,*) 'Integrando...'

      DO istep = 1, nsteps
!        --- Actualizar posiciones ---
!        x(t+dt) = x(t) + v(t) * dt
         DO i = 1, N
            x(i) = x(i) + vx(i) * dt
            y(i) = y(i) + vy(i) * dt
            z(i) = z(i) + vz(i) * dt
         END DO

!        --- Aplicar Condiciones de Contorno Periodicas ---
         DO i = 1, N
            IF (x(i) > L) x(i) = x(i) - L
            IF (x(i) < 0.0) x(i) = x(i) + L

            IF (y(i) > L) y(i) = y(i) - L
            IF (y(i) < 0.0) y(i) = y(i) + L

            IF (z(i) > L) z(i) = z(i) - L
            IF (z(i) < 0.0) z(i) = z(i) + L
         END DO

!        Las velocidades NO cambian (a = 0 para gas ideal)
      END DO

      WRITE(*,*) 'Integracion completada:', nsteps, 'pasos'

!======================================================================
!  SECCION 3.3: VERIFICACION
!======================================================================
!  Verificamos que las cantidades fisicas se conservan.
!  Para un gas ideal, energia y momento deben conservarse exactamente.
!======================================================================

!  --- Momento total ---
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))

      WRITE(*,*)
      WRITE(*,*) '=== VERIFICACION FINAL ==='
      WRITE(*,*) 'Momento total (debe ser ~0):'
      WRITE(*,*) 'Px =', Px
      WRITE(*,*) 'Py =', Py
      WRITE(*,*) 'Pz =', Pz

!  --- Energia ---
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      E_cinetica = 0.5 * suma_v2
      E_pot = 0.0
      E_total = E_cinetica + E_pot

      WRITE(*,*)
      WRITE(*,*) 'Energia cinetica final:', E_cinetica
      WRITE(*,*) 'Energia total final:', E_total
      WRITE(*,*) 'Diferencia energia:', ABS(E_total - 0.5 * suma_v2)

      WRITE(*,*)
      WRITE(*,*) 'NOTA: La energia NO ha cambiado porque'
      WRITE(*,*) 'las velocidades son constantes (gas ideal).'

!  --- Guardar datos ---
      OPEN(10, FILE='data/velocidades_dm.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) vx(i), vy(i), vz(i)
      END DO
      CLOSE(10)

      OPEN(10, FILE='data/modulos_dm.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) SQRT(vx(i)**2 + vy(i)**2 + vz(i)**2)
      END DO
      CLOSE(10)

      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - data/velocidades_dm.dat'
      WRITE(*,*) '  - data/modulos_dm.dat'

      END PROGRAM
