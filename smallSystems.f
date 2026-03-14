      PROGRAM smallSystems
      IMPLICIT NONE

!======================================================================
!  PROGRAMA v1: COLECTIVIDAD MICROCANONICA
!  Simulacion de N particulas con colisiones elasticas
!  Condiciones de contorno periodicas (PBC)
!  Autor: 
!  Fecha: 
!======================================================================

      INTEGER, PARAMETER :: NMAX = 2000
      INTEGER :: N, i, j, icoll, ncolisiones
      INTEGER :: i1, i2

!---------------------------------------------------------------------
!  Variables cinematicas
!---------------------------------------------------------------------
      REAL :: vx(NMAX), vy(NMAX), vz(NMAX)
      REAL :: x(NMAX), y(NMAX), z(NMAX)
      REAL :: v_modulo, vcm_x, vcm_y, vcm_z
      REAL :: L

!---------------------------------------------------------------------
!  Variables para angulos y numeros aleatorios
!---------------------------------------------------------------------
      REAL :: theta, phi, aleatorio

!---------------------------------------------------------------------
!  Variables termodinamicas
!---------------------------------------------------------------------
      REAL :: T_deseada, T_actual, suma_v2
      REAL :: Px, Py, Pz, E_total, E_inicial

!---------------------------------------------------------------------
!  Variables para colisiones
!---------------------------------------------------------------------
      REAL :: nx, ny, nz
      REAL :: v1x, v1y, v1z, v2x, v2y, v2z
      REAL :: prod_escalar
      REAL :: dx, dy, dz, rij

!======================================================================
!  SECCION 2.1: ESTRUCTURA BASICA - INICIALIZACION
!======================================================================

!  Parametros del sistema
      N = 1000
      v_modulo = 1.0
      ncolisiones = 200000
      L = 10.0

!---------------------------------------------------------------------
!  2.1.0: Inicializar posiciones aleatorias en el volumen
!---------------------------------------------------------------------
      DO i = 1, N
         x(i) = L * RAND()
         y(i) = L * RAND()
         z(i) = L * RAND()
      END DO

!---------------------------------------------------------------------
!  2.1.1: Inicializar velocidades con igual modulo y direccion aleatoria
!         Distribucion uniforme en esfera (angulos solidos iguales)
!---------------------------------------------------------------------
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
!  2.1.2: Imponer momento total cero (P = 0)
!         Restar velocidad del centro de masas a cada particula
!---------------------------------------------------------------------
      vcm_x = SUM(vx(1:N)) / REAL(N)
      vcm_y = SUM(vy(1:N)) / REAL(N)
      vcm_z = SUM(vz(1:N)) / REAL(N)

      DO i = 1, N
         vx(i) = vx(i) - vcm_x
         vy(i) = vy(i) - vcm_y
         vz(i) = vz(i) - vcm_z
      END DO

!======================================================================
!  SECCION 2.2: CONTROL TERMICO - ESCALADO DE VELOCIDADES
!======================================================================

      WRITE(*,*) '=== INICIALIZACION ==='
      WRITE(*,*) 'Verificacion momento total (inicial):'
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))
      WRITE(*,*) 'Px total =', Px
      WRITE(*,*) 'Py total =', Py
      WRITE(*,*) 'Pz total =', Pz

!---------------------------------------------------------------------
!  2.2.1: Calcular temperatura actual desde energia cinetica
!         T = <v^2> / 3  (unidades: m=kB=1)
!---------------------------------------------------------------------
      T_deseada = 1.5
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))

      WRITE(*,*)
      WRITE(*,*) 'Temperatura actual:', T_actual
      WRITE(*,*) 'Temperatura deseada:', T_deseada

!---------------------------------------------------------------------
!  2.2.2: Escalar velocidades para alcanzar temperatura deseada
!         v_i' = v_i * sqrt(T_deseada / T_actual)
!---------------------------------------------------------------------
      DO i = 1, N
         vx(i) = vx(i) * SQRT(T_deseada / T_actual)
         vy(i) = vy(i) * SQRT(T_deseada / T_actual)
         vz(i) = vz(i) * SQRT(T_deseada / T_actual)
      END DO

      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_inicial = 0.5 * suma_v2

      WRITE(*,*)
      WRITE(*,*) 'Temperatura despues de escalar:', T_actual
      WRITE(*,*) 'Energia cinetica inicial:', E_inicial

!======================================================================
!  SECCION 2.3: COLISIONES PARTICULA-PARTICULA
!======================================================================

!---------------------------------------------------------------------
!  2.3.1: Seleccionar pares de particulas aleatoriamente
!         Realizar ncolisiones colisiones elasticas
!         Colision elastica 3D con eje de colision aleatorio
!         Conserva momento total y energia cinetica
!---------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,*) '=== PARTICLE-PARTICLE COLLISIONS ==='

      DO icoll = 1, ncolisiones
!        Seleccionar dos particulas distintas
         i1 = INT(N * RAND()) + 1
         i2 = INT(N * RAND()) + 1
         DO WHILE (i2 == i1)
            i2 = INT(N * RAND()) + 1
         END DO

!        Generar eje de colision (vector unitario aleatorio)
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
         nx = SIN(theta) * COS(phi)
         ny = SIN(theta) * SIN(phi)
         nz = COS(theta)

!        Guardar velocidades antes de la colision
         v1x = vx(i1)
         v1y = vy(i1)
         v1z = vz(i1)
         v2x = vx(i2)
         v2y = vy(i2)
         v2z = vz(i2)

!        Formula colision elastica (masas iguales):
!        v1' = v1 - <v1-v2,n> * n
!        v2' = v2 + <v1-v2,n> * n
         prod_escalar = (v1x - v2x) * nx + (v1y - v2y) * ny + (v1z - v2z) * nz

         vx(i1) = v1x - prod_escalar * nx
         vy(i1) = v1y - prod_escalar * ny
         vz(i1) = v1z - prod_escalar * nz

         vx(i2) = v2x + prod_escalar * nx
         vy(i2) = v2y + prod_escalar * ny
         vz(i2) = v2z + prod_escalar * nz
      END DO

      WRITE(*,*) 'Colisiones completadas:', ncolisiones

!======================================================================
!  SECCION 2.4: CONDICIONES DE CONTORNO PERIODICAS (PBC)
!======================================================================

!---------------------------------------------------------------------
!  2.4.1: Actualizar posiciones con PBC
!         x = x + v*dt (simplificado: dt=1)
!         Si sale por un lado, entra por el opuesto
!         x' = x - L si x > L
!         x' = x + L si x < 0
!---------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,*) '=== PERIODIC BOUNDARY CONDITIONS ==='

      DO i = 1, N
         x(i) = x(i) + vx(i)
         y(i) = y(i) + vy(i)
         z(i) = z(i) + vz(i)

!        Aplicar PBC en cada direccion
         IF (x(i) > L) x(i) = x(i) - L
         IF (x(i) < 0.0) x(i) = x(i) + L

         IF (y(i) > L) y(i) = y(i) - L
         IF (y(i) < 0.0) y(i) = y(i) + L

         IF (z(i) > L) z(i) = z(i) - L
         IF (z(i) < 0.0) z(i) = z(i) + L
      END DO

      WRITE(*,*) 'Posiciones actualizadas con PBC'
      WRITE(*,*) 'Tamano de caja: L =', L

!======================================================================
!  SECCION 2.5: MUESTREO Y ANALISIS
!======================================================================

!---------------------------------------------------------------------
!  2.5.1: Verificar conservacion de magnitudes fisicas
!---------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,*) '=== VERIFICACION FINAL ==='

!        Momento total (debe conservarse exactamente con PBC)
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))
      WRITE(*,*) 'Momento total final:'
      WRITE(*,*) 'Px =', Px, '  Py =', Py, '  Pz =', Pz

!        Energia y temperatura (se conserva en colisiones elasticas)
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_total = 0.5 * suma_v2

      WRITE(*,*)
      WRITE(*,*) 'Temperatura final:', T_actual
      WRITE(*,*) 'Energia cinetica final:', E_total
      WRITE(*,*) 'Diferencia energia:', ABS(E_total - E_inicial)

!---------------------------------------------------------------------
!  2.5.2: Guardar datos para analisis (histograma, comparacion MB)
!---------------------------------------------------------------------
      OPEN(10, FILE='velocidades_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) vx(i), vy(i), vz(i)
      END DO
      CLOSE(10)

      OPEN(10, FILE='modulos_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) SQRT(vx(i)**2 + vy(i)**2 + vz(i)**2)
      END DO
      CLOSE(10)

      OPEN(10, FILE='posiciones_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) x(i), y(i), z(i)
      END DO
      CLOSE(10)

      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - velocidades_final.dat (vx, vy, vz por linea)'
      WRITE(*,*) '  - modulos_final.dat (|v| por linea)'
      WRITE(*,*) '  - posiciones_final.dat (x, y, z por linea)'

      END PROGRAM
