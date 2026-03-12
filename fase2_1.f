      PROGRAM inicializacion
      IMPLICIT NONE
      INTEGER, PARAMETER :: NMAX = 32
      INTEGER :: N, i
      REAL :: vx(NMAX), vy(NMAX), vz(NMAX)
      REAL :: v_modulo, vcm_x, vcm_y, vcm_z
      REAL :: theta, phi, aleatorio

      N = 20
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

      vcm_x = SUM(vx(1:N)) / REAL(N)
      vcm_y = SUM(vy(1:N)) / REAL(N)
      vcm_z = SUM(vz(1:N)) / REAL(N)

      DO i = 1, N
         vx(i) = vx(i) - vcm_x
         vy(i) = vy(i) - vcm_y
         vz(i) = vz(i) - vcm_z
      END DO

      WRITE(*,*) 'Velocidades después de imponer P=0:'
      DO i = 1, N
         WRITE(*,100) i, vx(i), vy(i), vz(i)
      END DO
 100  FORMAT('Particula ',I2,': vx=',F10.6,' vy=',F10.6,' vz=',F10.6)

      WRITE(*,*)
      WRITE(*,*) 'Verificacion momento total:'
      WRITE(*,*) 'Px total =', SUM(vx(1:N))
      WRITE(*,*) 'Py total =', SUM(vy(1:N))
      WRITE(*,*) 'Pz total =', SUM(vz(1:N))

      END PROGRAM
