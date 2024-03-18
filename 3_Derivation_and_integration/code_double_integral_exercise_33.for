program dobleintegral
  implicit none

  real*8 a, b, IT, x, alpha, ITreal, error
  integer i
  dimension alpha(0:100), x(0:100)

  print*, 'SHEET 3: DIFFERENTIATION AND INTEGRATION'
  print*, '3.3 - Surface Integral'
  print*

  print*, 'Enter a and b, the integration parameters in x:'
  read*, a, b
  print*, 'The integration parameters in y are: y1=0 and y2=1+x^2'

  ! Value obtained with Mathematica
  ITreal = 0.4747798814549d0

  ! Read the file with Gauss-Legendre weights
  open(200, file='cerosGaussLegendre.txt')
  do i = 0, 9
    read(200, *) x(i), alpha(i)
  end do

  call integral(x, alpha, a, b, IT)

  ! Calculate the percentage relative error
  error = ABS((IT - ITreal) / ITreal) * 1.d2

  print*, 'The value of the integral is:', IT
  print*, 'The percentage relative error is:', error
  print*, 'They have been saved in the file integraldoble.txt'

  ! Write results to file
  open(20, file='integraldoble.txt')
  write(20, '(10f15.10)') IT, error
  close(20)
  close(200)

  stop
end program

! Gauss-Legendre integration with ten points for x (integral in y already done)
subroutine integral(x, alpha, a, b, IT)
  implicit none
  real*8 z, IT, a, b, x, alpha, IY, aa, bb
  integer i
  dimension alpha(0:100), x(0:100)

  IT = 0.d0

  do i = 0, 9
    z = (b + a + (b - a) * x(i)) / 2.d0
    call integralY(aa, bb, x, alpha, z, IY)
    IT = IT + (b - a) * alpha(i) * IY / 2.d0
  end do

  return
end subroutine

! Gauss-Legendre integration with ten points for y
subroutine integralY(aa, bb, x, alpha, z, IY)
  implicit none
  real*8 z, IY, aa, bb, f, x, alpha, y
  integer j
  dimension alpha(0:100), x(0:100)

  aa = 0.d0
  bb = 1.d0 + z**2
  IY = 0.d0
  do j = 0, 9
    y = (bb + aa + (bb - aa) * x(j)) / 2.d0
    IY = IY + (bb - aa) * alpha(j) * f(z, y) / 2.d0
  end do

  return
end subroutine

! Function of two variables (x, y) to be integrated
function f(z, y)
  implicit none
  real*8 f, y, z

  f = (z**2 + y**2) * dexp(-z * y) * dcos(z) * dsin(y)

  return
end
