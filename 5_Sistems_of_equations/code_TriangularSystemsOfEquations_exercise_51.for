program sist_ec_triangulares
  implicit none
  real*8 :: matriz1(1:100,1:100), matriz2(1:100,1:100), b1(1:100), b2(1:100)
  integer :: i, j, n, tipo_matriz
  
  print*, 'SHEET 5: LINEAR SYSTEMS SOLUTION'
  print*, '5.1 - Triangular Systems of Equations'
  print*

  print*, 'Enter the dimension of the matrix:'
  read*, n

  open(10, file='matriz1.txt')
  open(20, file='matriz2.txt')
  open(15, file='b1.txt')
  open(16, file='b2.txt')

  if (n == 5) then
    do i = 1, n
      read(10,*) (matriz1(i,j), j = 1, n)
      read(15,*) b1(i)
    end do
  end if
  if (n == 6) then
    do i = 1, n
      read(20,*) (matriz2(i,j), j = 1, n)
      read(16,*) b2(i)
    end do
  end if

  print*, 'Is the matrix upper or lower triangular?'
  print*, 'Press [1] for upper; [2] for lower'

  read*, tipo_matriz
  print*, 'The solution is:'
  select case (tipo_matriz)
    case (1)
      call back_substitution(matriz2, b2, n)
    case (2)
      call forward_substitution(matriz1, b1, n)
    case default
      print*, 'Invalid case'
  end select

  stop
end program

subroutine back_substitution(matriz, b, n)
  implicit none
  real*8 :: x(1:100), b(1:100), matriz(1:100,1:100), suma, detA
  integer :: i, j, n

  do i = n, 1, -1
    suma = 0.d0
    do j = i+1, n
      suma = suma + matriz(i,j) * x(j)
    end do
    x(i) = (b(i) - suma) / matriz(i,i)
  end do

  do i = 1, n
    print*, x(i)
  end do

  detA = 1
  do i = n, 1, -1
    detA = detA * matriz(i,i)
  end do

  print*, 'The determinant is', detA
end subroutine

subroutine forward_substitution(matriz, b, n)
  implicit none
  real*8 :: x(1:100), b(1:100), matriz(1:100,1:100), suma, detA
  integer :: i, j, n

  do i = 1, n
    x(i) = 0.d0
  end do

  do i = 1, n
    suma = 0.d0
    do j = 1, i-1
      suma = suma + matriz(i,j) * x(j)
    end do
    x(i) = (b(i) - suma) / matriz(i,i)
    print*, x(i)
  end do

  detA = 1
  do i = 1, n
    detA = detA * matriz(i,i)
  end do

  print*, 'The determinant is', detA
end subroutine

