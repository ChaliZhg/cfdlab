      real*8 function decas(degree, coeff, t)
      implicit none
c
c    uses  de Casteljau to compute one coordinate
c    value of a  Bezier curve. Has to be called 
c    for each coordinate  (x,y, and/or z) of a control polygon.
c    Input:   degree: degree of curve.
c             coeff:  array with coefficients of curve.
c             t:      parameter value.
c     Output: coordinate value.
c
        integer degree
        real*8 coeff(*), t

        integer r,i
        real*8 t1
        real*8 coeffa(0:30)  ! an auxiliary array. Change dim. if too small

        t1 = 1.0d0 - t
        do i=0,degree  
            coeffa(i)=coeff(i+1) ! save input
        enddo

        do r=1,degree
           do i=0,degree - r
              coeffa(i)= t1* coeffa(i)  +   t * coeffa(i+1)
           enddo
        enddo

        decas = coeffa(0)
        end
