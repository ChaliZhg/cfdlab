C------------------------------------------------------------------------------
C.....Update the solution
C.....Current implementation is single step forward Euler
C------------------------------------------------------------------------------
      subroutine update(irk, divf, prim_old, prim, cvarea, spts, dt)
      implicit none
      include 'param.h'
      integer          irk, spts(nspmax)
      double precision divf(nvar,npmax), prim_old(nvar,npmax), 
     &                 prim(nvar,npmax), cvarea(npmax), dt(npmax)

      integer          i, j
      double precision con0_old(nvar), con1_old(nvar),
     &                 con_new(nvar,npmax)

      do i=1,np
         call prim_to_con(i, prim_old, con0_old)
         call prim_to_con(i, prim,     con1_old)
         do j=1,nvar-1
            con_new(j,i) = airk(irk)*con0_old(j) + 
     &                     birk(irk)*(con1_old(j) - 
     &                     (dt(i)/cvarea(i))*divf(j,i))
         enddo
      enddo
                  
      if(iflow .ne. inviscid)then
         do i=1,nsp
            j = spts(i)
            con_new(2,j) = 0.0d0
            con_new(3,j) = 0.0d0
         enddo
      endif

      call con_to_prim(con_new, prim)

      return
      end
