program polynomial_100

  implicit none
  real*8 x,f,sum,product,y,d
  dimension x(0:100),f(0:100)
  integer i,j,k,n,nodes
  n=18
  
  print*, 'SHEET 2: INTERPOLATION AND APPROXIMATION OF FUNCTIONS'
  print*, '2.2 - Lagrange vs. Splines (modified 2.1.b).'
  print*

  print*, 'Enter the number of nodes in the mesh (maximum 100):'
  read*, nodes
  ! Open the file with the given data
  open(33,file='data.txt')
  do i=0,n
    read(33,*) x(i),f(i)
  end do 
  close(33)

  d=(x(n)-x(0))/nodes ! Distance between the nodes
  do k=0,nodes
    ! Calculation of the distance between nodes
    y=x(0)+d*k 
    sum=0.d0 ! Initialize the summation

    do i=0,n
      product=1.d0 ! Initialize the product for each Li
      do j=0,n
        ! Calculate Li
        if (i.NE.j) then
          product=product*(y-x(j))/(x(i)-x(j))
        end if
      end do
      sum=sum+f(i)*product ! Lagrange polynomial
    end do
    ! Save the values of 'y' and 'sum' obtained
    open(44,file='lagrangeoutput.txt')
    write(44,*) y,sum
    close(44)
  end do

  stop 
end
