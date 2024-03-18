program lagrange

  implicit none
  real*8 x,f,sum,product,point
  dimension x(0:20),f(0:20)
  integer i,j
  integer n

  print*, 'SHEET 2: INTERPOLATION AND APPROXIMATION OF FUNCTIONS'
  print*, '2.1 - Lagrange Method. Part (a)'
  print*
  
  n=7
  ! Open the file with the given data
  open(34,file='data.txt')
  do i=0,n
    read(34,*) x(i),f(i)
  end do 
  close(34)

  print*,'Enter x* to calculate interpolated f(x*)'
  read*, point

  sum=0.d0 ! Initialize the summation

  do i=0,n
    product=1.d0 ! Initialize the product for each Li
    do j=0,n
      ! Calculate Li ensuring that x(i) is different from x(j)
      if (i.NE.j) then
        product=product*(point-x(j))/(x(i)-x(j))
      end if
    end do
    ! Sum all f(i)*Li to obtain the polynomial
    sum=sum+f(i)*product
  end do

  print*, 'The value of f(x*) is', sum

  stop 
end
