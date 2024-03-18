program bisection

  implicit none
  real*8 a,b,beta,d,c,old_c,f
  integer counter

  print*, 'SHEET 1: NONLINEAR UNIVARIABLE EQUATIONS'
  print*, '1.1 - Bisection Method.'
  print*

  ! Requesting problem data from the user
  print*,'Enter the interval.'
  read*, a,b
  d=f(a)*f(b)
  do while (d.GT.0.d0)
    print*,'Incorrect interval. Enter the interval.'
    read*, a,b
    d=f(a)*f(b)
  end do

  print*, 'Enter the precision (for convergence).'
  read*, beta

  ! Implementing the algorithm of the method
  c=(a+b)/2.d0
  old_c=0.d0
  counter=0
  do while (beta.LT.ABS(c-old_c))
    old_c=c
    d=f(a)*f(c)
    if (d.GT.0.d0) then
      a=c
    else
      b=c
    end if
    c=(a+b)/2.d0
    counter=counter+1
  end do

  ! Displaying the result
  print*, 'The solution is:', c
  print*, 'Number of iterations:', counter

  stop 
end

! Defining the function for which roots are to be found
function f(x)
  implicit none
  real*8 f,x
  
  f=0.5d0*dexp(-x/4.d0)*dsin(2.d0-(2.d0/5.d0)*x**2)
  
  return 
end
