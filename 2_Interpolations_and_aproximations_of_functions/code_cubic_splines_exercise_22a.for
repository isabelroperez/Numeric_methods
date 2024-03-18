program splines
      
  implicit none
  real*8 x,f,a,b,c,d,sx,ai,bi,ci,di,xi,h,y
  integer i,j,k,n,nodes
  dimension x(0:4000),f(0:4000),a(0:4000),b(0:4000),c(0:4000),
 & d(0:4000)

  print*, 'SHEET 2: INTERPOLATION AND APPROXIMATION OF FUNCTIONS'
  print*, '2.2 - Cubic Splines.'
  print*

  ! Passing the data from the .txt file to the program
  open(10, file='data.txt')
  open(20, file='output.txt')
  
  write(*,*) 'Enter the number of input data points'
  read(*,*) n
  n=n-1 ! Counting starts at 0, so the number of data points is n-1
  write(*,*) 'Enter the number of mesh points (max 100)'
  read(*,*) nodes
  
  do i=0,n
    read(10,*) x(i),f(i) ! Reading from the input file
  end do
     
  h=(x(n)-x(0))/nodes  ! Interval distance
  k=0 ! Initialize k
  
  ! Call the subroutine to calculate the coefficients
  call coef(x,f,n,a,b,c,d)
  
  do j=0,nodes
    ! Calculating the position of y between two x points
    y=x(0)+h*j
    do i=0,n
      if ((y.GT.x(i)).AND.(y.LT.x(i+1))) then
        k=i
      end if        
    end do      
    ! After finding the coefficients, construct the function by joining the splines
    ai=a(k)
    bi=b(k)
    ci=c(k)
    di=d(k)
    xi=x(k)
    sx=ai*(y-xi)**3+bi*(y-xi)**2+ci*(y-xi)+di ! Interpolation at each point
    write(20,*) y,sx
  end do

  close(10)
  close(20)
  write(*,*) 'Interpolation data saved successfully.'      

  stop
end

! SUBROUTINE TO CALCULATE THE COEFFICIENTS
subroutine coef(xx,ff,nn,aa,bb,cc,dd)
  implicit none
  real*8 xx,ff,aa,bb,cc,dd,A,B,C,r,hh
  integer i,nn
  dimension xx(0:4000),ff(0:4000),aa(0:4000),bb(0:4000),cc(0:4000),
 & dd(0:4000),A(0:4000),B(0:4000),C(0:4000),r(0:4000),hh(0:4000)

  do i=0,nn-1
    hh(i)=xx(i+1)-xx(i)
  end do
  do i=0,nn
    dd(i)=ff(i)
  end do
  do i=1,nn-1
    r(i)=(3.d0/hh(i))*(dd(i+1)-dd(i))-(3.d0/hh(i-1))*(dd(i)-dd(i-1))
    B(i)=2.d0*(hh(i)+hh(i-1))
  end do
  do i=2,nn-1
    A(i)=hh(i-1)
  end do
  do i=1,nn-2
    C(i)=hh(i)
  end do
  do i=2,nn-1
    B(i)=B(i)-(C(i-1)*A(i))/B(i-1)
    r(i)=r(i)-(r(i-1)*A(i))/B(i-1)
  end do
  bb(nn-1)=r(nn-1)/B(nn-1)
  do i=nn-2,1,-1
    bb(i)=(r(i)-C(i)*bb(i+1))/B(i)
  end do
  bb(0)=0
  bb(nn)=0
  do i=0,nn-1
    cc(i)=(dd(i+1)-dd(i))/(xx(i+1)-xx(i))-
 & ((bb(i+1)+2.d0*bb(i))/3.d0)*(xx(i+1)-xx(i))
  end do
  do i=1,nn
    aa(i-1)=(bb(i)-bb(i-1))/(3.d0*(xx(i)-xx(i-1)))
  end do 

  return
end subroutine
