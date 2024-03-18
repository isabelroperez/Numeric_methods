program LU_decomposition

  implicit none
  real*8 :: a(1:100,1:100), b3(1:100), m(1:100,1:100), u(1:100,1:100), y(1:100), x(1:100), suma, detA
  integer :: i, j, n, p, k

  print*, 'SHEET 5: LINEAR SYSTEMS SOLUTION'
  print*, '5.2 - LU Decomposition for Non-Triangular Systems'
  print*

  print*, 'Enter the dimension of the matrix:'
  read*, n

  open(10, file='a.txt')
  open(16, file='b3.txt')    
  open(11, file='m.txt')
  open(12, file='u.txt')
  open(3, file='detA.txt')

  do i = 1, n
    do j = 1, n
      a(i, j) = 0.d0
      m(i, j) = 0.d0
      u(i, j) = 0.d0
    end do
    b3(i) = 0.d0
  end do

  do i = 1, n
    read(10,*) (a(i, j), j = 1, n)
    read(16,*) b3(i)
  end do

  do i = 1, n
    m(i, i) = 1.d0
    u(1, i) = a(1, i)
  end do

  do k = 2, n
    do i = 1, k - 1
      suma = 0.d0
      do j = 1, k - 1
        suma = suma + m(k, j) * u(j, i)
      end do

      if (ABS(u(i, i)).GT.1.D-10) then
        m(k, i) = (a(k, i) - suma) / u(i, i)            
      else 
        print*, 'Error in denominator'
        stop
      end if
    end do

    do p = k, n
      suma = 0.d0
      do j = 1, n - 1
        suma = suma + m(k, j) * u(j, p)
      end do
      u(k, p) = a(k, p) - suma
    end do
  end do

  do i = 1, n
    write(11, '(10f9.4)') (m(i, j), j = 1, n)
    write(12, '(10f9.4)') (u(i, j), j = 1, n)
  end do

  detA = 1
  do i = 1, n
    detA = detA * u(i, i)
  end do

  print*, 'The determinant is', detA
  Write(3,*) detA

  call forward_substitution(m, b3, n, y)
  call backward_substitution(u, y, n, x)

  close(3)
  close(11)
  close(12)
  close(10)
  close(16)

  stop
end program

subroutine forward_substitution(m, b3, n, y)
  implicit none
  real*8 :: y(1:100), b3(1:100), m(1:100,1:100), suma
  integer :: i, j, n

  do i = 1, n
    y(i) = 0.d0
  end do

  do i = 1, n
    suma = 0.d0
    do j = 1, i - 1
      suma = suma + m(i, j) * y(j)
    end do
    y(i) = (b3(i) - suma) / m(i, i)
  end do
end subroutine

subroutine backward_substitution(u, y, n, x)
  implicit none
  real*8 :: x(1:100), y(1:100), u(1:100,1:100), suma
  integer :: i, j, n

  do i = n, 1, -1
    suma = 0.d0
    do j = i + 1, n
      suma = suma + u(i, j) * x(j)
    end do
    x(i) = (y(i) - suma) / u(i, i)
  end do

  print*, 'The resulting vector is:'
  do i = 1, n
    print*, x(i)
  end do
end subroutine
