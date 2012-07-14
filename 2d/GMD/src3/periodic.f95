subroutine periodic(con)
   use comvar
   implicit none

   real :: con(nvar,-1:nx+2,-1:ny+2,-1:nz+2)

   integer :: i, j, k

   do j=1,ny
      do k=1,nz
         con(:,-1,  j,k) = con(:,nx-1,j,k)
         con(:, 0,  j,k) = con(:,nx,  j,k)
         con(:,nx+1,j,k) = con(:,1,   j,k)
         con(:,nx+2,j,k) = con(:,2,   j,k)
      enddo
   enddo

   do i=1,nx
      do k=1,nz
         con(:,i,  -1,k) = con(:,i,ny-1,k)
         con(:,i,   0,k) = con(:,i,ny  ,k)
         con(:,i,ny+1,k) = con(:,i,1   ,k)
         con(:,i,ny+2,k) = con(:,i,2   ,k)
      enddo
   enddo

   do i=1,nx
      do j=1,ny
         con(:,i,j,  -1) = con(:,i,j,nz-1)
         con(:,i,j,   0) = con(:,i,j,nz  )
         con(:,i,j,nz+1) = con(:,i,j,1   )
         con(:,i,j,nz+2) = con(:,i,j,2   )
      enddo
   enddo

   ! Corners cells

   ! Near k=0
   do i=-1,0
      do j=-1,0
         do k=-1,0
            con(:,i,j,k) = con(:,i+nx,j+ny,k+nz)
         enddo
      enddo
   enddo

   do i=nx+1,nx+2
      do j=-1,0
         do k=-1,0
            con(:,i,j,k) = con(:,i-nx,j+ny,k+nz)
         enddo
      enddo
   enddo

   do i=nx+1,nx+2
      do j=ny+1,ny+2
         do k=-1,0
            con(:,i,j,k) = con(:,i-nx,j-ny,k+nz)
         enddo
      enddo
   enddo

   do i=-1,0
      do j=ny+1,ny+2
         do k=-1,0
            con(:,i,j,k) = con(:,i+nx,j-ny,k+nz)
         enddo
      enddo
   enddo

   ! Near k=nz
   do i=-1,0
      do j=-1,0
         do k=nz+1,nz+2
            con(:,i,j,k) = con(:,i+nx,j+ny,k-nz)
         enddo
      enddo
   enddo

   do i=nx+1,nx+2
      do j=-1,0
         do k=nz+1,nz+2
            con(:,i,j,k) = con(:,i-nx,j+ny,k-nz)
         enddo
      enddo
   enddo

   do i=nx+1,nx+2
      do j=ny+1,ny+2
         do k=nz+1,nz+2
            con(:,i,j,k) = con(:,i-nx,j-ny,k-nz)
         enddo
      enddo
   enddo

   do i=-1,0
      do j=ny+1,ny+2
         do k=nz+1,nz+2
            con(:,i,j,k) = con(:,i+nx,j-ny,k-nz)
         enddo
      enddo
   enddo

end subroutine periodic
