C.....Read parameters from an input file and set freestream values
      subroutine read_input
      implicit none
      include 'param.h'
      integer inp, iargc, n, inpstatus
      character sdummy*32

      n = iargc()
      if(n .eq. 0)then
            print*,'You must specify an input file.'
            stop
      endif

      call getarg(1,inpfile)

      inp = 11
      open(unit=inp, file=inpfile, status='old')
      print*,'Reading parameters from ',inpfile
      read(inp,*)sdummy, istart
      read(inp,*)sdummy, iflow
      read(inp,*)sdummy, mach_inf
      read(inp,*)sdummy, aoa_deg
      read(inp,*)sdummy, Rey
      read(inp,*)sdummy, cfl
      read(inp,*)sdummy, iterlast
      read(inp,*)sdummy, maxiter
      read(inp,*)sdummy, minresidue
      read(inp,*)sdummy, saveinterval
      read(inp,*)sdummy, niso
      read(inp,*)sdummy, iflux
      read(inp,*)sdummy, ILIMIT
      read(inp,*)sdummy, vortex, xref, yref
      read(inp,*)sdummy, cell_type
      read(inp,*)sdummy, gridfile
      close(inp)

      inpstatus = yes

      if(istart .ne. scratch .and. istart .ne. restart)then
            print*,'Unknown start option',istart
            print*,'Possible values: 1=scratch or 2=restart'
            inpstatus = no
      endif

      if(iflow .ne. inviscid .and. iflow .ne. laminar .and.
     &   iflow .ne. turbulent)then
            print*,'Unknown flow type',iflow
            print*,'Possible values: 1=inviscid, 2=laminar, 3=turbulent'
            inpstatus = no
      endif

      if(iflux .ne. iroe .and. iflux .ne. ikfvs .and. 
     &   iflux .ne. ihllc)then
            print*,'Unknown flux',iflux
            print*,'Possible values: 1=roe, 2=kfvs'
            inpstatus = no
      endif

      if(ilimit .ne. no .and. ilimit .ne. yes)then
            print*,'Unknown limiter option',ilimit
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(vortex .ne. yes .and. vortex .ne. no)then
            print*,'Unknown vortex option',vortex
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(cell_type .ne. median .and. cell_type .ne. barth)then
            print*,'Unknown cell type',cell_type
            print*,'Possible values: 1=median cell, 2=barth cell'
            inpstatus = no
      endif

      if(inpstatus .eq. no) stop

      return
      end
