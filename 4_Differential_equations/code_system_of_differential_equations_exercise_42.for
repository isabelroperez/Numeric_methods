program sist_ec_dif
  implicit none
  real*8 :: h, vi, yi, t, datos, v0, t0, y0, k1, k2
  integer :: n, ii, i
  character(len=50) :: a, b

  print*, 'SHEET 4: ORDINARY DIFFERENTIAL EQUATIONS'
  print*, '3.2 - Vertical Projectile Motion (System of Differential Equations)'
  print*

  print*, 'Use predefined values or input your own?'
  print*, 'Press [1] for predefined values:'
  read*, datos

  if (datos.EQ.1) then
    h = 1.d-3
    t = 1.6
    v0 = 8.d0
    t0 = 0.d0
    y0 = 0.d0
    k2 = 1.d-2
    k1 = 4.d-3   
  else
    print*, 'Please enter the flight time:'
    read*, t

    print*, 'Please enter the step size:'
    read*, h
  
    print*, 'Please enter the initial velocity:'
    read*, v0

    print*, 'Please enter the initial height:'
    read*, y0

    print*, 'Finally, enter the initial time:'
    read*, t0

    print*, 'Enter the values of k1 and k2:'
    read*, k1, k2
  end if 

  n = t / h

  print*, 'Calculating data for two different k values'

  do i = 0, 1
    select case (i)
      case (0)
        a = 'altitude_k1.txt'
        b = 'velocity_k1.txt'
        ii = 100
        call rk4(h, n, vi, yi, v0, t0, y0, k1, a, b, ii)
      case (1)
        a = 'altitude_k2.txt' 
        b = 'velocity_k2.txt'
        ii = 102
        call rk4(h, n, vi, yi, v0, t0, y0, k2, a, b, ii)
    end select
  end do
  
  stop
end

subroutine rk4(h, n, vi, yi, v0, t0, y0, k, a, b, ii)
  implicit none
  real*8 :: yi, yvieja, ymax, t0, ti, tmax, tchoque, v0, vi, h, Ky1, Ky2, Ky3, Ky4, Kv1, Kv2, Kv3, Kv4, k
  integer :: i, n, ii
  character(len=50) :: a, b

  open(ii, file=a)
  open(ii+1, file=b)

  do i = 0, n-1
    ti = t0 + i * h
    yvieja = yi

    Ky1 = fy(ti, yi, vi)
    Kv1 = fv(ti, yi, vi, k) 
    Ky2 = fy(ti+h/2.d0, yi+h/2.d0*Ky1, vi+h/2.d0*Kv1)
    Kv2 = fv(ti+h/2.d0, yi+h/2.d0*Ky1, vi+h/2.d0*Kv1, k)
    Ky3 = fy(ti+h/2.d0, yi+h/2.d0*Ky2, vi+h/2.d0*Kv2)
    Kv3 = fv(ti+h/2.d0, yi+h/2.d0*Ky2, vi+h/2.d0*Kv2, k)
    Ky4 = fy(ti+h, yi+h*Ky3, vi+h*Kv3)
    Kv4 = fv(ti+h, yi+h*Ky3, vi+h*Kv3, k)
        
    yi = yi + h/6.d0 * (Ky1 + 2.d0*Ky2 + 2.d0*Ky3 + Ky4)
    vi = vi + h/6.d0 * (Kv1 + 2.d0*Kv2 + 2.d0*Kv3 + Kv4)
    write(ii, *) ti, yi
    write(ii+1, *) ti, vi

    if (yvieja < yi) then
      ymax = yi
      tmax = ti
    end if

    if (yi > 0.d0) then
      tchoque = ti
    end if
  end do
      
  print*, 'Maximum altitude and time reached:', ymax, tmax
  print*, 'Flight time:', tchoque
  close(ii)
  close(ii+1)
end subroutine

function fy(t, y, v)
  implicit none
  real*8 :: t, y, v, fy

  fy = v

  return
end function

function fv(t, y, v, k)
  implicit none
  real*8 :: t, y, v, fv, g, m, fy
  real*8, parameter :: g = 9.807d0, m = 0.11d0

  fy = v
  fv = -k/m * fy(t, y, v) * abs(fy(t, y, v)) - g

  return
end function
