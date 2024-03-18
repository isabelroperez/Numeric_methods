program newton

  implicit none
  real*8 g,h,a,b,beta,bold
  integer c,nmax

  print*, 'SHEET 1: NONLINEAR UNIVARIABLE EQUATIONS'
  print*, '1.3 - Newton's Method.'
  print*,

  print*, 'Enter a starting point.'
  read*, a
  do while (ABS(h(a)).LT.ABS(1.d-16))
    print*, 'Incorrect point. Enter a starting point.'
    read*, a
  end do
  
  ! Using a maximum number of iterations to avoid infinite loop         
  print*, 'Enter a maximum number of iterations.'
  read*, nmax
  
  print*, 'Enter the precision (for convergence)'
  read*, beta

  ! Starting the algorithm of the method
  b=a-g(a)/h(a)
  bold=0.d0
  c=0
  do while ((beta.LT.ABS(b-bold)) ! Convergence condition
 & .AND. (ABS(h(a)).GT.ABS(1.d-16)).AND.nmax.GT.c) ! Division by zero and infinite loop condition
    bold=b
    a=b
    b=a-g(a)/h(a)
    c=c+1
  end do
  
  if (nmax.NE.c) then
    print*, 'The root found is', b
    print*, 'Number of iterations:', c
  else
    print*, 'No solution found.'
  end if
  
  stop
end

! Working with the first and second derivative; we need the relative extremes of the function
function g(x) ! First derivative
  implicit none
  real*8 g,x
  
  g=(-1.d0/40.d0)*dexp(-x/4.d0)*(16.d0*x*dcos(2.d0-(2.d0/5.d0)*x**2)
 &+5.d0*dsin(2.d0-(2.d0/5.d0)*x**2))
  
  return 
end

function h(x) ! Second derivative
  implicit none
  real*8 h,x
  
  h=(1.d0/800.d0)*dexp(-x/4.d0)*(160.d0*(x-2.d0)*dcos(2.d0
 &-(2.d0/5.d0)*x**2)+(25.d0-256.d0*x**2)*dsin(2.d0
 &-(2.d0/5.d0)*x**2))
  
  return 
end
