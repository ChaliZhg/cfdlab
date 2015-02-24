subroutine save_along_line(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, q, a, m
   character(len=512) :: filename

   if(test_case.eq.ibriowu)then
      j = ny/2
   else if(test_case.eq.ialfven)then
      j = 1
   else
      return
   endif

   filename = 'line.dat'
   open(10,file=trim(filename))

   do i=1,nx
      x = xmin + (i-1)*dx + 0.5*dx

      q = sqrt(pri(2,i,j)**2 + pri(3,i,j)**2 + pri(4,i,j)**2)
      a = sqrt(gamma*pri(5,i,j)/pri(1,i,j))
      m = q/a
      write(10,'(10E24.14)') x, pri(1,i,j), pri(2,i,j), pri(3,i,j), &
                             pri(4,i,j), pri(5,i,j), pri(6,i,j), &
                             pri(7,i,j), pri(8,i,j), m
   enddo

   close(10)

end subroutine save_along_line
