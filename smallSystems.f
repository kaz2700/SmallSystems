      PROGRAM smallSystems
      IMPLICIT NONE

      INTEGER, PARAMETER :: NMAX = 2000
      INTEGER :: N, i, icoll, ncolisiones
      INTEGER :: i1, i2

!  DECLARACION DE VARIABLES
      
!  vx, vy, vz: componentes de velocidad de cada particula
!  x, y, z:     posiciones de cada particula
      REAL :: vx(NMAX), vy(NMAX), vz(NMAX)
      REAL :: x(NMAX), y(NMAX), z(NMAX)

!  Parametros del sistema
      REAL :: v_modulo         ! Velocidad inicial de cada particula
      REAL :: vcm_x, vcm_y, vcm_z  ! Velocidad del centro de masas
      REAL :: L               ! Lado del cubo (condiciones periodicas)

!  Variables para numeros aleatorios
      REAL :: theta, phi, aleatorio

!  Variables termodinamicas
      REAL :: T_deseada       ! Temperatura objetivo
      REAL :: T_actual        ! Temperatura actual (calculada)
      REAL :: suma_v2         ! Sumatoria de v^2
      REAL :: Px, Py, Pz      ! Componentes del momento total
      REAL :: E_total          ! Energia total
      REAL :: E_inicial        ! Energia inicial (para verificar conservacion)

!  Variables para colisiones elasticas
      REAL :: nx, ny, nz      ! Componentes del eje de colision (vector unitario)
      REAL :: v1x, v1y, v1z   ! Velocidad de particula 1 antes de colision
      REAL :: v2x, v2y, v2z   ! Velocidad de particula 2 antes de colision
      REAL :: prod_escalar     ! Producto escalar (v1-v2)·n

!  SECCION 2.1

      N = 1000                ! Numero de particulas
      v_modulo = 1.0         ! Velocidad inicial (magnitud)
      ncolisiones = 10   ! Numero de colisiones a simular
      L = 10.0               ! Tamano de la caja

      DO i = 1, N
         x(i) = L * RAND()
         y(i) = L * RAND()
         z(i) = L * RAND()
      END DO

!  Para una distribucion uniforme sobre la esfera unitaria, usamos:
!  - theta = arccos(2*u - 1) con u~U(0,1)  -> distribucion uniforme en [0, pi]
!  - phi = 2*pi*u con u~U(0,1)             -> distribucion uniforme en [0, 2pi]
      DO i = 1, N
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio

         vx(i) = v_modulo * SIN(theta) * COS(phi)
         vy(i) = v_modulo * SIN(theta) * SIN(phi)
         vz(i) = v_modulo * COS(theta)
      END DO

!  En la colectividad microcanonica, el momento total debe ser cero.
!  v_i' = v_i - V_cm
!  donde V_cm = (1/N) * sum(v_i)
      vcm_x = SUM(vx(1:N)) / REAL(N)
      vcm_y = SUM(vy(1:N)) / REAL(N)
      vcm_z = SUM(vz(1:N)) / REAL(N)

      DO i = 1, N
         vx(i) = vx(i) - vcm_x
         vy(i) = vy(i) - vcm_y
         vz(i) = vz(i) - vcm_z
      END DO

!  SECCION 2.2
!  La temperatura se relaciona con la energia cinetica media:
!    T = <v^2> / 3   (en unidades donde m = kB = 1)
!
!  Para cambiar la temperatura de T1 a T2:
!    v' = v * sqrt(T2/T1)

      WRITE(*,*) 'INICIALIZACION'
      WRITE(*,*) 'Verificacion momento total (inicial):'

!  Calcular momento total inicial
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))
      WRITE(*,*) 'Px total =', Px
      WRITE(*,*) 'Py total =', Py
      WRITE(*,*) 'Pz total =', Pz

      T_deseada = 1.5
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 /  REAL(N) / 3.0

      WRITE(*,*)
      WRITE(*,*) 'Temperatura actual (antes de escalar):', T_actual
      WRITE(*,*) 'Temperatura deseada:', T_deseada

      DO i = 1, N
         vx(i) = vx(i) * SQRT(T_deseada / T_actual)
         vy(i) = vy(i) * SQRT(T_deseada / T_actual)
         vz(i) = vz(i) * SQRT(T_deseada / T_actual)
      END DO

!  Verificar que se alcanza la temperatura deseada
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_inicial = 0.5 * suma_v2

      WRITE(*,*)
      WRITE(*,*) 'Temperatura despues de escalar:', T_actual
      WRITE(*,*) 'Energia cinetica inicial:', E_inicial

!  SECCION 2.3
!  Despues de suficientes colisiones, el sistema alcanza equilibrio
!  termico y la distribucion de velocidades es Maxwell-Boltzmann.

!  Algoritmo de colision elastica para masas iguales:
!  1. Seleccionar dos particulas al azar
!  2. Generar un eje de colision aleatorio (vector unitario n)
!  3. Calcular el producto escalar (v1-v2)·n
!  4. Actualizar velocidades:
!       v1' = v1 - ((v1-v2)·n) * n
!       v2' = v2 + ((v1-v2)·n) * n

      WRITE(*,*)
      WRITE(*,*) 'SIMULANDO COLISIONES'
      WRITE(*,*) 'Numero de colisiones:', ncolisiones

      DO icoll = 1, ncolisiones

!  Seleccionar dos particulas distintas

         i1 = INT(N * RAND()) + 1
         i2 = INT(N * RAND()) + 1
         DO WHILE (i2 == i1)
            i2 = INT(N * RAND()) + 1
         END DO

!        --- Generar eje de colision aleatorio ---
!        Usamos el mismo metodo que para velocidad inicial
         aleatorio = RAND()
         theta = ACOS(2.0 * aleatorio - 1.0)
         aleatorio = RAND()
         phi = 2.0 * 3.141592653589793 * aleatorio
         nx = SIN(theta) * COS(phi)
         ny = SIN(theta) * SIN(phi)
         nz = COS(theta)

!        --- Guardar velocidades antes de la colision ---
         v1x = vx(i1)
         v1y = vy(i1)
         v1z = vz(i1)
         v2x = vx(i2)
         v2y = vy(i2)
         v2z = vz(i2)

!        --- Aplicar formula de colision elastica ---
         prod_escalar = (v1x - v2x) * nx + (v1y - v2y) * ny + (v1z - v2z) * nz

         vx(i1) = v1x - prod_escalar * nx
         vy(i1) = v1y - prod_escalar * ny
         vz(i1) = v1z - prod_escalar * nz

         vx(i2) = v2x + prod_escalar * nx
         vy(i2) = v2y + prod_escalar * ny
         vz(i2) = v2z + prod_escalar * nz
      END DO

      WRITE(*,*) 'Colisiones completadas'

!  SECCION 2.4

!  Las condiciones de contorno periodicas simulan un sistema infinito. Si una particula cruza un borde, aparece por el opuesto.

      WRITE(*,*)
      WRITE(*,*) 'APLICANDO CONDICIONES DE CONTORNO PERIODICAS'

!  Actualizar posiciones: x = x + v * dt (aquí dt = 1)
      DO i = 1, N
         x(i) = x(i) + vx(i)
         y(i) = y(i) + vy(i)
         z(i) = z(i) + vz(i)
      END DO

!  Aplicar condiciones en cada direccion
      DO i = 1, N
         IF (x(i) > L) x(i) = x(i) - L
         IF (x(i) < 0.0) x(i) = x(i) + L

         IF (y(i) > L) y(i) = y(i) - L
         IF (y(i) < 0.0) y(i) = y(i) + L

         IF (z(i) > L) z(i) = z(i) - L
         IF (z(i) < 0.0) z(i) = z(i) + L
      END DO

      WRITE(*,*) 'Posiciones actualizadas con condiciones de contorno periodicas'

!  SECCION 2.5: VERIFICACION Y GUARDADO DE DATOS

!  Verificamos que las cantidades fisicas se conservan correctamente.
!  Guardamos los datos para analisis posterior (histogramas, etc.)
!======================================================================

      WRITE(*,*)
      WRITE(*,*) '=== VERIFICACION FINAL ==='

!  --- Verificar conservacion del momento total ---
!  Debe ser aproximadamente cero (error de precision numerica)
      Px = SUM(vx(1:N))
      Py = SUM(vy(1:N))
      Pz = SUM(vz(1:N))
      WRITE(*,*) 'Momento total final:'
      WRITE(*,*) 'Px =', Px, '  Py =', Py, '  Pz =', Pz

!  --- Verificar conservacion de energia ---
!  La energia debe mantenerse constante (a menos por errores numericos)
      suma_v2 = SUM(vx(1:N)**2 + vy(1:N)**2 + vz(1:N)**2)
      T_actual = suma_v2 / (3.0 * REAL(N))
      E_total = 0.5 * suma_v2

      WRITE(*,*)
      WRITE(*,*) 'Temperatura final:', T_actual
      WRITE(*,*) 'Energia cinetica final:', E_total
      WRITE(*,*) 'Diferencia energia (inicial - final):', ABS(E_total - E_inicial)

!  --- Guardar datos para analisis ---
      OPEN(10, FILE='data/velocidades_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) vx(i), vy(i), vz(i)
      END DO
      CLOSE(10)

      OPEN(10, FILE='data/modulos_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) SQRT(vx(i)**2 + vy(i)**2 + vz(i)**2)
      END DO
      CLOSE(10)

      OPEN(10, FILE='data/posiciones_final.dat', STATUS='REPLACE')
      DO i = 1, N
         WRITE(10,*) x(i), y(i), z(i)
      END DO
      CLOSE(10)

      WRITE(*,*)
      WRITE(*,*) 'Archivos de datos generados:'
      WRITE(*,*) '  - data/velocidades_final.dat (vx, vy, vz por linea)'
      WRITE(*,*) '  - data/modulos_final.dat (|v| por linea)'
      WRITE(*,*) '  - data/posiciones_final.dat (x, y, z por linea)'

      END PROGRAM
