program derivative

  implicit none
  real*8 f, x, h, h0, deriv1, deriv2, rich, value, e1, e2, e3
  integer i, n, j
  character(len=50) :: method, ext, filename

  print*, 'SHEET 3: DIFFERENTIATION AND INTEGRATION'
  print*, '3.1 - Differentiation'
  print*
  
  ! Saving data in two files, one for each value of y(x)
  do j = 0, 1
    select case (j)
      ! Values in each case are calculated with Mathematica
      case (0)
        print*
        print*, 'Case 1: y(x) = x^3'
        value = 0.492431480686d0
        n = 1
      case (1)
        print*
        print*, 'Case 2: y(x) = 1 + x^2'
        value = -0.849593314370d0
        n = 2
    end select
    
    method = 'derivatives'
    ext = '.txt'
    write(filename, '(A9,I1,A4)') method, n, ext
    open(50, file=filename)

    ! Initial values
    x = 1.d0
    h0 = 0.01d0
    print*, 'Results are displayed in this order:'
    print*, 'Non-symmetric, symmetric, Richardson extrapolation,'
    print*, 'followed by their respective errors:'

    do i = 0, 9
      ! Calculate h using i (without subtracting from the previous one)
      h = h0 - i * 1.d-3
      
      ! Non-symmetric formula (definition)
      deriv1 = (y(j, x + h) - y(j, x)) / h
      
      ! Symmetric formula
      deriv2 = (y(j, x + h) - y(j, x - h)) / (2.d0 * h)
      
      ! Richardson extrapolation
      rich = -1.d0 / (6.d0 * h) * (y(j, x + h) - 8.d0 * y(j, x + h / 2.d0) -
     & y(j, x - h) + 8.d0 * y(j, x - h / 2.d0))
      
      ! Calculate relative errors of each calculation in percentage
      e1 = ABS((deriv1 - value) / value) * 1.d2
      e2 = ABS((deriv2 - value) / value) * 1.d2
      e3 = ABS((rich - value) / value) * 1.d2

      ! Specify the precision of each number
      write(*, 300) deriv1, deriv2, rich, e1, e2, e3
300   format(3(f12.8,1x), 3(e16.8,1x))

      write(50, '(10f15.10)') deriv1, deriv2, rich, e1, e2, e3
    end do
    close(50)
  end do

  stop
end

! Define functions
function y(j, x)
  implicit none
  real*8 y, x
  integer j     
  
  select case (j)
    case (0)
      y = x**3
    case (1)
      y = 1 + x**2
  end select

  return 
end function y

function f(j, x)
  implicit none
  real*8 f, x
  integer j

  f = dexp(-y(j, x) / 1.d1) * dcos((dlog(y(j, x)) + 3.d0))**2
  
  return 
end function f
