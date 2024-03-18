program secant

  implicit none
  real*8 f,a,b,d,beta,c,old_c
  integer counter,nmax

  print*, 'SHEET 1: NONLINEAR UNIVARIABLE EQUATIONS'
  print*, '1.2 - Secant Method.'
  print*,

  print*, 'Enter two points for the secant.'
  read*, a,b
  d=f(a)-f(b)
  do while (ABS(d).LT.ABS(1.d-15)) ! In this case, the method diverges (division by 0)
    print*, 'Invalid values. Enter two points for the secant.'
    read*, a,b
  end do

  ! Introducing nmax to avoid infinite loop
  print*, 'Enter a maximum number of iterations.'
  read*, nmax     
  print*, 'Enter the precision (for convergence).'
  read*, beta

  ! Execution of the method's algorithm
  c=b-f(b)*(b-a)/(f(b)-f(a))
  old_c=0.d0
  counter=0
  do while ((beta.LT.ABS(c-old_c)) ! Convergence condition
 & .AND. (ABS(f(b)-f(a)).GT.ABS(1.d-16)).AND.nmax.GT.counter) ! Infinite loop and division by zero condition
    old_c=c
    a=b
    b=c
    c=b-f(b)*(b-a)/(f(b)-f(a))
    counter=counter+1
  end do

  ! Displaying results
  if (nmax.NE.counter) then
    print*, 'The root found is', c
    print*, 'Number of iterations:', counter
  else
    print*, 'No solution found.'
  end if
  
  stop
end

! Definition of the function
function f(x)
  implicit none
  real*8 f,x
  
  f=0.5d0*dexp(-x/4.d0)*dsin(2.d0-(2.d0/5.d0)*x**2)
  
  return 
end

