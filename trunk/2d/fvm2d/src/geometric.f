C Calculate geometric quantities like face areas, normals, volumes, etc.
      subroutine geometric(ptype, elem, edge, coord, elarea, cvarea, ds,
     &                     dsb, cvareamc, wd, drmin, spts, fpts, opts,
     &                     bpts)
      implicit none
      include 'param.h'

      integer          ptype(npmax), elem(nvemax,ntmax),
     &                 esup1(mesup*npmax), esup2(npmax+1),
     &                 psup1(mpsup*npmax), psup2(npmax+1),
     &                 edge(2,nemax), edneigh(2,nemax), spts(nspmax),
     &                 fpts(nfpmax), opts(nopmax), bpts(nbpmax)
      double precision coord(2, npmax), elarea(ntmax), cvarea(npmax), 
     &                 ds(2,nemax), dsb(2,npmax), drmin(npmax), 
     &                 tcoord(2,ntmax), cvareamc(npmax), wd(npmax)

c     Make sure elements are oriented ccw
      call tri_orient(elem, coord)

c     Find elements surrounding a point
      call el_surr_point(elem, esup1, esup2)

c     Find points surrounding a point
      call pt_surr_pt(esup1, esup2, elem, psup1, psup2)

c     Create edges
      call create_edge(psup1, psup2, edge)

c     Find element adjoining each edge
      call el_surr_edge(esup1, esup2, elem, edge, edneigh)

c     Write grid in gnuplot format for visualization
      call write_grid(coord, edge, edneigh)

c     Calculate element and control volume areas
      call areas(coord, tcoord, elem, elarea, cvarea, cvareamc)

c     Find wall distance for Spalart-Allmaras model
      if(iflow .eq. turbulent) call wall_dist(edge, ptype, coord, spts,
     &                                        psup1, psup2, wd)

c     Find length of control volume faces
      call flength(coord, tcoord, edge, edneigh, ds, dsb)

c     Length scale for time-step calculation
      call dtlength(coord, elarea, elem, drmin)

      return
      end
