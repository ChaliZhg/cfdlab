
      allocate(mmat(nc),amat(nc),bmat(nc),cmat(nc))

       do i=1,nc
         dt(i) = CFL * dx(i)/1.5
         mmat(i) = dx(i)/dt(i)
       enddo

      dxl = xc(1) - xv(1)
      dxr = xc(2) - xc(1)
      bmat(1) = -(1.0/dxl + 1.0/dxr)
      cmat(1) =  1.0/dxr
      do i=2,nc-1
         dxl = xc(i) - xc(i-1)
         dxr = xc(i+1) - xc(i)
         amat(i) =  1.0/dxl
         bmat(i) = -(1.0/dxl + 1.0/dxr)
         cmat(i) =  1.0/dxr
      enddo
      dxl = xc(nc) - xc(nc-1)
      dxr = xv(nc+1) - xc(nc)
      amat(nc) =  1.0/dxl
      bmat(nc) = -(1.0/dxl + 1.0/dxr)

      amat(1:nc) = -amat(1:nc)
      bmat(1:nc) = mmat(1:nc) - bmat(1:nc)
      cmat(1:nc) = -cmat(1:nc)
