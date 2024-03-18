program integral

  implicit none
  real*8 a, b, G, GL, s, TNA, TA, SNA, r, es, eG, eGL, eTNA, eTA, eSNA

  print*, 'SHEET 3: DIFFERENTIATION AND INTEGRATION'
  print*, '3.2 - Simple Integral'
  print*

  ! Calculate the exact value with Mathematica
  r = 0.6092797476926595d0
  
  print*, 'Enter a and b (integration parameters):'
  read*, a, b

  print*, 'Numerical integration calculation for each method:'
  print*

  ! Format specification
  300   format(f12.8,1x,e16.8)

  ! Simpson's Rule
  call simpson(a, b, r, s, es)
  print*, 'Simpson and its error:'
  write(*, 300) s, es
  print*, '--------------------------------'

  ! Newton-Cotes Closed with 6 equispaced points (n=5)
  call newtoncotes(a, b, r, G, eG)
  print*, 'Newton-Cotes and its error:'
  write(*, 300) G, eG
  print*, '--------------------------------'

  ! Gauss-Legendre with 10 points (i=0,9, n=9)
  call gausslegendre(a, b, r, GL, eGL)
  print*, 'Gauss-Legendre and its error:'
  write(*, 300) GL, eGL
  print*, '--------------------------------'

  ! Non-adaptive Trapezoidal with 32 points (n=31)
  call trapecioNA(a, b, r, TNA, eTNA)
  print*, 'Non-adaptive Trapezoidal and its error:'
  write(*, 300) TNA, eTNA
  print*, '--------------------------------'

  ! Adaptive Trapezoidal for convergence beta=0.00001d0
  call trapecioA(a, b, r, TA, eTA)
  print*, 'Adaptive Trapezoidal and its error:'
  write(*, 300) TA, eTA
  print*, '--------------------------------'

  ! Non-adaptive Simpson with 51 points (2m=50)
  call simpsonNA(a, b, r, SNA, eSNA)
  print*, 'Non-adaptive Simpson and its error:'
  write(*, 300) SNA, eSNA
  print*, '--------------------------------'

  ! Write results to file
  open(30, file='integrals.txt')
  write(30, 300) s, es
  write(30, 300) G, eG
  write(30, 300) GL, eGL
  write(30, 300) TNA, eTNA
  write(30, 300) TA, eTA
  write(30, 300) SNA, eSNA
  close(30)

  stop
end

! Simpson's Rule
subroutine simpson(aa, bb, r, s, es)
  implicit none
  real*8 aa, bb, medio, h, f, s, r, es
  medio = (aa + bb) / 2.d0
  h = (bb - aa) / 2.d0
  s = h / 3.d0 * (f(aa) + 4.d0 * f(medio) + f(bb))
  es = ABS((s - r) / r) * 1.d2
  return
end

! Newton-Cotes Closed
subroutine newtoncotes(aa, bb, r, G, eG)
  implicit none
  real*8 A, aa, bb, d, Gvieja, G, f, h, r, eG
  dimension A(0:100)
  integer i, n
  n = 5
  A(0) = 0.0989583d0
  A(1) = 0.390625d0
  A(2) = 0.260417d0
  A(3) = 0.260417d0
  A(4) = 0.390625d0
  A(5) = 0.0989583d0
  h = (bb - aa) / n
  Gvieja = 0
  do i = 0, n
    d = aa + i * h
    G = A(i) * f(d)
    G = G + Gvieja
    Gvieja = G
  end do
  eG = ABS((G - r) / r) * 1.d2
  return
end

! Gauss-Legendre
subroutine gausslegendre(aa, bb, r, GL, eGL)
  implicit none
  real*8 aa, bb, z, x, alpha, GL, f, r
