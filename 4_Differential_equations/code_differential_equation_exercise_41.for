program ec_diferencial
  implicit none
  real*8 :: x0, xn, y0, h
  integer :: i, n

  print*, 'SHEET 4: ORDINARY DIFFERENTIAL EQUATIONS'
  print*, '3.1 - Differential Equation'
  print*

  print*, 'Calculating solution to differential equation.'
  print*, 'No action required, please wait.'
  print*, 'Files will be saved in the current directory.'
  print*, 'Generated files:'

  ! Initial condition and mesh (same for all methods)
  y0 = 0.d0
  x0 = 0.d0
  xn = 4.d0

  ! Automate the change of mesh
  do i = 0, 3
    select case (i)
      case (0)
        n = 10
      case (1)
        n = 25
      case (2)
        n = 50
      case (3)
        n = 100
      case default
        print*, 'Error occurred with mesh selection.'
    end select
    h = (xn - x0) / (n - 1)
    call taylor(x0, y0, h, n)
    call heun(x0, y0, h, n)
    call oilar(x0, y0, h, n)
    call rk4(x0, y0, h, n)
  end do

  print*;
  print*, 'Successfully completed.'
  
  stop
end program

! Differential equation function
function f(x, y)
  implicit none
  real*8 :: x, y, f

  f = y - dexp(x) * dcos(x) * dsin(x)

  return
end function

! Taylor method
subroutine taylor(x0, y0, h, n)
  implicit none
  real*8 :: x0, xi, y0, yi, h
  integer :: i, n
  character(len=50) :: method, ext, filename

  method = 'taylor'
  ext = '.txt'

  do i = 0, n - 1
    xi = x0 + i * h
    yi = y0
    open(50, file=trim(method)//trim(adjustl(str(i,3)))//trim(ext))
    do
      write(50, *) xi, yi
      if (xi >= x0 + (n - 1) * h) exit
      yi = yi + h * f(xi, yi)
      xi = xi + h
    end do
    close(50)
  end do

  return
end subroutine

! Heun's method
subroutine heun(x0, y0, h, n)
  implicit none
  real*8 :: x0, xi, y0, yi, h, xnueva
  integer :: i, n
  character(len=50) :: method, ext, filename

  method = 'heun'
  ext = '.txt'

  do i = 0, n - 1
    xi = x0 + i * h
    yi = y0
    open(60, file=trim(method)//trim(adjustl(str(i,3)))//trim(ext))
    do
      write(60, *) xi, yi
      xnueva = x0 + (i + 1) * h
      if (xi >= x0 + (n - 1) * h) exit
      yi = yi + h / 2.d0 * (f(xi, yi) + f(xnueva, yi + h * f(xi, yi)))
      xi = xi + h
    end do
    close(60)
  end do

  return
end subroutine

! Modified Euler's method (Midpoint method)
subroutine oilar(x0, y0, h, n)
  implicit none
  real*8 :: x0, xi, y0, yi, h
  integer :: i, n
  character(len=50) :: method, ext, filename

  method = 'oilar'
  ext = '.txt'

  do i = 0, n - 1
    xi = x0 + i * h
    yi = y0
    open(70, file=trim(method)//trim(adjustl(str(i,3)))//trim(ext))
    do
      write(70, *) xi, yi
      if (xi >= x0 + (n - 1) * h) exit
      yi = yi + h * f(xi + h / 2.d0, yi + h / 2.d0 * f(xi, yi))
      xi = xi + h
    end do
    close(70)
  end do

  return
end subroutine

! Runge-Kutta 4th order method
subroutine rk4(x0, y0, h, n)
  implicit none
  real*8 :: x0, xi, y0, yi, h, K1, K2, K3, K4
  integer :: i, n
  character(len=50) :: method, ext, filename

  method = 'rk4'
  ext = '.txt'

  do i = 0, n - 1
    xi = x0 + i * h
    yi = y0
    open(80, file=trim(method)//trim(adjustl(str(i,3)))//trim(ext))
    do
      write(80, *) xi, yi
      if (xi >= x0 + (n - 1) * h) exit
      K1 = f(xi, yi)
      K2 = f(xi + h / 2.d0, yi + h / 2.d0 * K1)
      K3 = f(xi + h / 2.d0, yi + h / 2.d0 * K2)
      K4 = f(xi + h, yi + h * K3)
      yi = yi + h / 6.d0 * (K1 + 2.d0 * K2 + 2.d0 * K3 + K4)
      xi = xi + h
    end do
    close(80)
  end do

  return
end subroutine
