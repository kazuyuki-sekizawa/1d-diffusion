!--------+---------+---------+---------+---------+---------+---------+--
!        10        20        30        40        50        60        70
!--------+---------+---------+---------+---------+---------+---------+--
!! @file 1d-diffuion.f90
!> @brief Solves a 1D diffusion equation with the FTCS method.
!> @author Kazuyuki Sekizawa
!> @date 2022/05/20
!> @see https://nuclphystitech.wordpress.com/undergrad/comptphys/step-3/
!>
!> @details
!>   This program solves a one-dimensional diffusion equation,
!>   \f[
!>      \frac{\mathrm{d}f(x,t)}{\mathrm{d}t} = \kappa\frac{\mathrm{d}^2f(x,t)}{\mathrm{d}x^2}
!>   \f]
!>   using a forward difference for time derivative and a central difference
!>   for spatial derivative.
!>
!> @param i           Spatial index
!> @param ip1         Spatial index for i+1
!> @param im1         Spatial index for i-1
!> @param iter        Time index
!> @param Niter_max   The maximum number of time steps to be computed
!> @param dx          Grid spacing
!> @param dt          Time step
!> @param kappa       Diffusion coefficient
!> @param f(1:Nx)     An array of the function to be evolved in time
!> @param f_new(1:Nx) A buffer for the array f
!> @param char_iter   Stores 'iter' as a 4-digit character, e.g., 0001, 0002, ...
!>
!> @pre Modify the input file, "input.txt".
!> @pre Specify the initial condition for f(x)
!>
!> <b>Compile</b>
!> @code $ gfortran 1d-diffusion.f90
!> @endcode
!>
!> <b>Execution</b>
!> @code $ ./a.out
!> @endcode
!>
!> @todo Test various doxygen functions
!> @todo Try Github
!>
!> @note    You can add note here!
!>
!> @warning You can add warning here!
!--------+---------+---------+---------+---------+---------+---------+--

PROGRAM diffusion_1d
  IMPLICIT NONE
  INTEGER :: i,ip1,im1,iter,Niter_max,Nx
  DOUBLE PRECISION :: dx!> dx
  DOUBLE PRECISION :: dt!< dt
  DOUBLE PRECISION :: kappa
  DOUBLE PRECISION,ALLOCATABLE :: f(:),f_new(:)
  CHARACTER(4) :: char_iter

  ! Read parameters
  OPEN(100,file="input.txt",status="old")
  READ(100,*) Nx,dx
  READ(100,*) Niter_max,dt
  READ(100,*) kappa

  ALLOCATE(f(Nx),f_new(Nx))

  ! Initial condition
  f=0d0; f(Nx/2)=10d0

  ! Output the initial condition
  OPEN(100,file="f_0000.dat")
  DO i=1,Nx
    WRITE(100,*) i*dx, f(i)
  END DO
  CLOSE(100)

  ! Loop over time
  DO iter=1,Niter_max
    ! Time-evolution
    DO i=1,Nx
      ip1=i+1; im1=i-1
      IF(i==1 ) im1=Nx
      IF(i==Nx) ip1=1
      f_new(i)=f(i)+kappa*dt/dx**2*(f(ip1)-2d0*f(i)+f(im1))
    END DO
    f=f_new

    ! Output
    WRITE(char_iter,"(i0.4)") iter
    OPEN(100,file="f_"//char_iter//".dat")
    DO i=1,Nx
      WRITE(100,*) i*dx, f(i)
    END DO
    CLOSE(100)
    !...
  END DO

  DEALLOCATE(f,f_new)
  STOP
END PROGRAM diffusion_1d

!--------+---------+---------+---------+---------+---------+---------+--
