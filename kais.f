!#####################################################################
! This program calculate the kai parameters for poor-solvent 
! using MC integration for a 2D nanopore
!
!#####################################################################

      subroutine kai

      use mvariables
      use mparameters
      use mkai
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mkaps

      implicit none

      include 'mpif.h'
      include 'MPI.h'

      real*8  suma(dimR, dimR+1, -Xulimit:Xulimit)  ! Xu(a,b, c)
                                                     ! interaction of a
                                                     ! segment at R = a
                                                     ! and z = 0, with a
                                                     ! segment at R = b
                                                     ! and Z = c
      real*8 R,theta,z
      real*8 thetamin, thetamax, rmax, rmin, zmax, zmin, r0
      real*8 cutoff
      real*8 cutoffmin

      real*8 rn
      integer i, ii
      integer iR, iZ, itheta
      real*8 radio ! total radius
      real*8 D, V ! total radius
      real*8 x1,x2,y1, y2, z1, z2, vect
      integer Rj, Zj, j
      character*40 filename
      real*8 sumXuii
      real*8 rounding

      radio = delta*float(Kapd)
      D = 2.0*radio
      if(readkai.ne.1) then
         if(rank.eq.0)print*,'Kai calculation, readkai =', readkai
         write(filename,'(A5, I4.4, A4)')
     & 'kais-', dimR, '.kai'
         open(unit=111, file=filename)
         suma = 0.0
         Xu = 0.0               ! vector Xu
         seed = readseed

         cutoff = float(Xulimit)*delta
         cutoffmin = 2.0*radio+delta

         rounding = 1e-6 ! used to prevent errors due to rounding

         do ii = 1, dimR       ! loop over segment positions

         zmin = -cutoff+rounding
         zmax = cutoff-rounding
            
         r0 = (dfloat(ii) - 0.5)*delta

         if(r0.gt.cutoff) then
         rmin = r0-cutoff+rounding
         rmax = r0+cutoff-rounding
         thetamax = asin(cutoff/r0)
         thetamin = -thetamax
         else
         rmin = 0+rounding
         rmax = r0+cutoff-rounding
         thetamin = 0
         thetamax = 2*pi
         endif


         do iR = 1, MCsteps
         do iZ = 1, MCsteps
         do itheta = 1, MCsteps
       
            Z = (zmax-zmin)*float(iZ-1)/float(MCsteps-1) + zmin
            theta = 
     & (thetamax-thetamin)*float(itheta-1)/float(MCsteps-1) + thetamin
            R = (rmax-rmin)*float(iR-1)/float(MCsteps-1) + rmin

!     segment coordinates (x1,y1,z1) and point to integrate coordinates
!     (x2,y2,z2)
               x1 = r0 ! theta segment = 0, z segment =
               y1 = 0.0
               z1 = 0.0
               x2 = R*cos(theta)
               y2 = R*sin(theta)
               z2 = z
               vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diff
               Rj = int(R/delta)+1 ! R cell
               Zj = nint(z/delta) ! from -.5 to .5 => layer 0
                                ! from .5 to 1.5 => layer 1
                                ! from -.5 to -1.5 => layer -1
                                ! z is the distance to the segment,
                                ! while R
                                ! is the position in the lattice
               if(Zj.gt.Xulimit) then
                print*, z, z/delta,i,ii
                stop
               endif

               if(Rj.gt.dimR)Rj=dimR+1 ! goes to bulk
               suma(ii, Rj, Zj) = suma(ii, Rj, Zj) + R

               if(vect.gt.cutoff)cycle ! outside cut-off sphere
               if(vect.lt.cutoffmin)cycle

               V = vect

               Xu(ii, Rj, Zj) = Xu(ii, Rj, Zj) 
     &        + R/12.0*((D*D)/((2*D+V)*V)+
     &        (D*D)/((D+V)*(D+V))+2.0*log((2*D+V)*V/((D+V)*(D+V))))

!               print*, ii, Rj, Zj, Xu(ii,Rj,Zj)

!               Xu(ii, Rj, Zj) = Xu(ii, Rj, Zj) + ((lseg/vect)**6)*R ! incluye el jacobiano R(segmento)

            enddo ! iR
            enddo ! iZ
            enddo ! itheta

            do Rj = 1, dimR+1
            do Zj = -Xulimit,Xulimit
              Xu(ii, Rj, Zj) = Xu(ii, Rj,Zj)/(MCsteps**3)*(zmax-zmin)*
     &                 (thetamax-thetamin)*(rmax-rmin)

              suma(ii, Rj, Zj) =suma(ii,Rj,Zj)/(MCsteps**3)*(zmax-zmin)*
     &                 (thetamax-thetamin)*(rmax-rmin)
            enddo
            enddo

            if(ii.eq.1) then
            sumXu = 0.0
            do Rj = 1,dimR+1
            do Zj = -Xulimit,Xulimit
            sumXu = sumXu+Xu(1,Rj,Zj) ! normalization for Xu
            write(111,*)ii,Rj,Zj,Xu(ii,Rj,Zj)
            enddo 
            enddo       
            else
            sumXuii = 0.0
            do Rj = 1,dimR+1
            do Zj = -Xulimit,Xulimit
            sumXuii = sumXuii+Xu(ii,Rj,Zj) 
            enddo 
            enddo
!            print*, ii, sumXu,sumXuii
            do Rj = 1,dimR+1
            do Zj = -Xulimit,Xulimit
            Xu(ii,Rj,Zj) = Xu(ii,Rj,Zj)*sumXu/sumXuii
            write(111,*)ii,Rj,Zj,Xu(ii,Rj,Zj)
            enddo
            enddo
            endif

         end do                 ! ii
         close(111)
      else                      ! readkai = 1

        write(filename,'(A11, I4.4, A4)')'../../kais-', dimR, '.kai'
         open(unit=111, file=filename)
         if(rank.eq.0)print*,'Read Kais'
        do j = 1, dimR*(dimR+1)*3
            read(111,*)ii,Rj,Zj,Xu(ii,Rj,Zj)
         enddo
      endif

! OJO 

! determine sum of kais for bulk

      sumXu = 0.0
      do Rj = 1,dimR+1
      do Zj = -Xulimit,Xulimit
      sumXu = sumXu+Xu(1,Rj,Zj) ! all should be the same, except the last ones.
      enddo 
      enddo       

      if(rank.eq.0)print*, 'sumXu', sumXu

      end

      subroutine lookup_kai

      use mlookup
      use mparameters
      use mkai
      use mncells
      use mporesystem
      use mparameters_chain
      use mparameters_monomer
        

      implicit none
      integer iC, iR, iZ,  iiR, iiZ
      integer iiiZ

! define Xu

      nXu = 0
      Xulist_cell = 0
      Xulist_value = 0

      do iC = 1, ncells ! loop over cells

      iR = indexa(iC,1) ! position of cell in the system
      iZ = indexa(iC,2)

      do iiR = 1, dimR+1 ! loop over R-neighbors
      do iiZ = -Xulimit,Xulimit ! loop over Z-neighbors

      if(Xu(iR,iiR,iiZ).ne.0d0) then ! Xu 

      iiiZ = iiZ+iZ

      do while (iiiZ.gt.dimZ)
      iiiZ = iiiZ-dimZ
      enddo
      
      do while (iiiZ.lt.1)
      iiiZ = iiiZ+dimZ
      enddo

      if (iiR.eq.dimR+1) then ! bulk
         nXu(iC) = nXu(iC) + 1
         if(nXu(iC).gt.size(Xulist_cell,2)) then
          print*, size(Xulist_cell,2)
          print*, 'Error in lookup_kai, increase Xulist_cell dimensions'
         stop
         endif
      Xulist_cell(iC, nXu(iC)) = ncells+1
      Xulist_value(iC, nXu(iC)) = Xu(iR, iiR, iiZ)

      else ! not bulk

        if(matriz(iiR, iiiZ).gt.0) then ! not a wall

        nXu(iC) = nXu(iC) + 1
        if(nXu(iC).gt.size(Xulist_cell,2)) then
         print*, size(Xulist_cell,2)
         print*, 'Error in lookup_kai, increase Xulist_cell dimensions'
         stop
        endif

      Xulist_cell(iC, nXu(iC)) = matriz(iiR, iiiZ)
      Xulist_value(iC, nXu(iC)) = Xu(iR, iiR, iiZ)

       endif ! not a wall
      endif ! bulk
      endif ! Xu =! 0

      enddo ! iiZ
      enddo ! iiR

      enddo ! iC
      return
      end



