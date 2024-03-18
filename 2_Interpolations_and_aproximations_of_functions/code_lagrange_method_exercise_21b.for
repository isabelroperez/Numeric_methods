program lagrange2

  implicit none
  real*8 x,f,sum,product,y,d
  dimension x(0:10),f(0:10)
  integer i,j,k,n,nodes
  n=7

  print*, 'SHEET 2: INTERPOLATION AND APPROXIMATION OF FUNCTIONS'
  print*, '2.1 - Lagrange Method. Part (b)'
  print*

  print*, 'Enter the number of nodes in the mesh (maximum 100):'
  read*, nodes

  ! Open the file with the given data
  open(34,file='data.txt')
  do i=0,n
    read (34,*) x(i),f(i)
  end do 
  close(34)

  open(45,file='meshdata.txt') ! Open the file for output data
  
  d=(x(n)-x(0))/nodes ! Distance between the nodes
  do k=0,nodes
    ! Added for part (b)
    y=x(0)+d*k ! Calculation of the distance between nodes
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
    write (45,*) y,sum
  end do

  close(45)
  
  stop 
end
