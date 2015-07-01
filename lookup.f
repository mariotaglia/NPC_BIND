      subroutine lookup

      use mlookup
      use mparameters
      use mporesystem
      use mncells
      use mparameters_chain
      use mparameters_monomer

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
      integer iR, iZ, indice, iC

c ----------------- contruye matriz -------------------------


C codes matriz = ! NCELLS+1: bulk iz
                 ! NCELLS+2: bulk derecho
                 ! -1:symmetry
                 ! -2: wall
                 !  0: canal
                 ! -3: no definido


c------------------ crea lookup list --------------------------

! indexa matriz


      indice = 0

      do iR = 1,dimR
      do iZ = 1,dimZ

      if(matriz(iR, iZ).eq.0) then
      indice = indice + 1
       matriz(iR, iZ) = indice
      indexa(indice,1) = iR
      indexa(indice,2) = iZ
      endif
      enddo
      enddo

! save to disk

c      do iR = 0,dimR+1
c      do iZ = 0,dimZ+1
c      write(998, *)iR, iZ, matriz(iR, iZ)
c      enddo
c      enddo

! Ahora define las matrices...


      do iR = 1,dimR
      do iZ = 1,dimZ

      if(matriz(iR,iZ).gt.0) then ! punto dentro del sistema

! Radio +1

      select case (matriz(iR+1,iZ))
 
      case (1 : )
      rp(matriz(iR,iZ))=matriz(iR+1,iZ)
      case (-2:-1)
      rp(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ
      
      end select

! Radio -1

      select case (matriz(iR-1,iZ))

      case (1 : )
      rm(matriz(iR,iZ))=matriz(iR-1,iZ)
      case (-2:-1)
      rm(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select

! Zeta +1

      select case (matriz(iR,iZ+1))

      case (1 : )
      zp(matriz(iR,iZ))=matriz(iR,iZ+1)
      if(iZ.eq.dimZ)zp(matriz(iR,iZ))=matriz(iR,1)
      case (-2:-1)
      zp(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select

! Zeta - 1

      select case (matriz(iR,iZ-1))

      case (1 : )
      zm(matriz(iR,iZ))=matriz(iR,iZ-1)
      if(iZ.eq.1)zm(matriz(iR,iZ))=matriz(iR,dimZ)
      case (-2:-1)
      zm(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select
      endif

      enddo
      enddo
  
      return
      end


       subroutine geom
       use mparameters
       use mncells
       use mporesystem
       use mparameters_chain
       use mparameters_monomer
       use mKaps

       implicit none
       integer midR, iR, iZ, midZ

       integer*4 temp(0:dimR+1, 0:dimZ+1)
       real*8 zpos, rpos, slope, ord, rcurve, rcurve2

       real*8 b_param, a_param, c_param
       real*8 volprot

       integer Zmin, Zmax

      integer PposZ

      common /PposZ/ PposZ
      


C codes matriz = ! NCELLS+1: bulk iz
                 ! NCELLS+2: bulk derecho
                 ! -1:symmetry
                 ! -2: wall
                 !  0: canal
                 ! -3: no definido


       matriz = -3 ! -3 => undefined

       ncells = 0
       temp = 0

!!!!!!!!!!! pore !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      select case (poretype)

! -----------------------  Simplest case: long pore --------------------------
      case (1)

       do iZ = 1, dimZ
       do iR = 1, CdimR
  
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1

       enddo
       enddo

! ---------------------------- empty ------------------------------------
      case (0)

       do iZ = 1, dimZ
       do iR = 1, dimR
  
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1

       enddo
       enddo


!------------------------ Hour glass, using a parabola ---------------------- 
       case (3)
       a_param = curvature
       b_param = -a_param*dimz*delta
       c_param = CdimR*delta
     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta

       do iZ = 1, dimZ
       do iR = 1, dimR

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

       rcurve = a_param*zpos**2 + b_param*zpos + c_param

       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ)).
     & or.((dimR-iR).le.Kapd)) 
     & then
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1
       endif

       enddo
       enddo

! ---------------------------------- conic pore ----------------------
       case (4)

       do iZ = 1, dimZ
       do iR = 1, dimR

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

       slope =
     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
       ord = dfloat(TdimR)*delta-RdimZ*delta*slope      

       rcurve = slope*zpos + ord

       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ))
     &.or.(iR.gt.dimR-2*Kapd-1)) then
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1
       endif

       enddo
       enddo


! ------------------------------ hour glass (parabola) and protein --------
c       case (5)
c
c       matriz = 0
c
c       temp = 0 
c       do iZ = 1, dimZ
c       temp(0,iZ) = -1
c       enddo
c       ncells = dimZ*dimR
c
c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       do iZ = 1, dimZ
c       do iR = 1, dimR
c
c       zpos = (dfloat(iZ)-0.5)*delta
c       rpos = (dfloat(iR)-0.5)*delta
c
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       if((rcurve.le.rpos).and.(iZ.gt.RdimZ).and.(iZ.le.(dimZ-RdimZ)))
c     & then ! es pared
c
c       matriz(iR, iZ) = -3 
c       ncells = ncells - 1
c
c      endif
c
c       if((zpos.gt.(dfloat(PposZ)-dfloat(PdimR))*delta)
c     & .and.(zpos.lt.(dfloat(PposZ)+dfloat(PdimR))*delta)) then
c
c       rcurve2 = sqrt((float(PdimR)*delta)**2
cc     &- (zpos-dfloat(PposZ)*delta)**2) 
c
c      if(rcurve2.gt.rpos) then ! es particula
c
c      if(matriz(iR, iZ).eq.-3) then
c      print*, 'Colision entre pared y particula'
c      endif
c
c      protein(iR,iZ) = prot_vol
c      temp(iR,iZ) = 1
c
c      endif
c      endif
c      
c      enddo
c      enddo
c
c      do iR = 1, dimR
c      do iZ = 1, dimZ
c
c      if(temp(iR,iZ).eq.1) then
c
c      if((temp(iR+1,iZ).eq.0).or.(temp(iR-1,iZ).eq.0).or.
c     & (temp(iR,iZ+1).eq.0).or.(temp(iR,iZ-1).eq.0)) then
c
c      volprot = volprot + delta**3*2*pi*(dfloat(iR)-0.5)
c      proteinq(iR,iZ) = prot_q
c
c      endif
c      endif
c     
c      enddo
c      enddo
c
c      if(volprot.ne.0)proteinq = proteinq / volprot
c
c       case (40)
c
c       do iZ = 1, dimZ
c       do iR = 1, dimR
c
c       zpos = (dfloat(iZ)-0.5)*delta
c       rpos = (dfloat(iR)-0.5)*delta
c
c       slope =
c     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
c       ord = dfloat(TdimR)*delta-RdimZ*delta*slope
c
c       rcurve = slope*zpos + ord
c
c       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ)))
c     & then
c       matriz(iR, iZ) = 0 ! canal
c       ncells = ncells + 1
c       endif
c
c      enddo
c      enddo

!------------------------ Hour glass, using half circle  ---------------------- 
       case (10)
     
       ncells = dimZ*dimR
       matriz = 0

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       do iZ = zmin+1, zmax
       do iR = 1, dimR

       rcurve = sqrt((dfloat(CdimZ)/2*delta)**2 
     & - ((dfloat(iZ-zmin)-0.5)*delta - dfloat(CdimZ)/2*delta)**2)
       rcurve = CdimR*delta + dfloat(CdimZ)/2*delta - rcurve 

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

!       if(rcurve.lt.rpos)
       if((rcurve.lt.rpos).and.((dimR-iR).gt.(2*Kapd)))
     & then
       matriz(iR, iZ) = -3 
       ncells = ncells - 1
       endif

       enddo
       enddo

!------------------------ Hour glass like  ---------------------- 

         case (11)
  
         ncells = dimZ*dimR
         matriz = 0
 
         zmin = (dimZ - CdimZ)/2
         zmax = zmin + CdimZ
 
         do iZ = zmin+1, zmin+CdimZ/2
         do iR = CdimR+1, dimR
 
         slope = (dfloat(CdimZ-CdimZmin)/2)/(dfloat(CdimZ)/2)
         slope = -1/slope
         ord = (dfloat(CdimR)+dfloat(CdimZ)/2)*delta-zmin*delta*slope
 
         zpos = (dfloat(iZ)-0.5)*delta
         rpos = (dfloat(iR)-0.5)*delta

         rcurve = slope*zpos + ord

         if((rcurve.lt.rpos).and.((dimR-iR).gt.2*Kapd))
     &   then
         matriz(iR, iZ) = -3
         ncells = ncells - 1
         endif
 
         enddo
         enddo

         do iZ = zmin+CdimZ/2+1, zmax
         do iR = CdimR+1, dimR

         slope = (dfloat(CdimZ-CdimZmin)/2)/(dfloat(CdimZ)/2)
         slope = 1/slope
         ord = (dfloat(CdimR)+dfloat(CdimZ)/2)*delta-zmax*delta*slope

         zpos = (dfloat(iZ)-0.5)*delta
         rpos = (dfloat(iR)-0.5)*delta
         rcurve = slope*zpos + ord
         if((rcurve.lt.rpos))
     &   then
         matriz(iR, iZ) = -3
         ncells = ncells - 1
         endif

         enddo
         enddo

!------------------------------ Short Cycilnder --------------------

         case (12)

         ncells = dimZ*dimR
         matriz = 0

         zmin = int((dimZ - CdimZ)/2)
         zmax = zmin + CdimZ

         do iZ = zmin+1, zmax
         do iR = CdimR+1, dimR-2*Kapd-1

         matriz(iR, iZ) = -3
         ncells = ncells - 1

         enddo
         enddo

       endselect


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reservoirs

       
       midZ = int(dfloat(dimZ)/2)
       do iZ = 0, midZ
       matriz(dimR+1, iZ) = ncells+1
       enddo

       do iZ = midZ, dimZ+1
       matriz(dimR+1, iZ) = ncells+2
       enddo

       do iR = 1, dimR+1
       matriz(iR, 0) = ncells+1
       matriz(iR, dimZ+1) = ncells+2
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! symmetry

       do iZ = 0, dimZ+1
       matriz(0, iZ) = -1
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! search for borders
       temp = matriz

       do iR = 1, dimR
       do iZ = 1, dimZ

          if(temp(iR, iZ).eq.-3) then
               if((temp(iR, iZ+1).eq.0).or.(temp(iR, iZ-1).eq.0).or.
     &   (temp(iR+1, iZ).eq.0).or.(temp(iR-1, iZ).eq.0)) then

                matriz(iR, iZ) = -2 ! pared
           endif
       endif

       enddo
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save to disk


c      do iR = 0,dimR+1
c      do iZ = 0,dimZ+1
c      write(999, *)iR, iZ, matriz(iR, iZ)
c      enddo
c      enddo
      end

     
      subroutine graftpoints
      use mparameters
      use mncells
      use mporesystem
      use mgraft       
      use mparameters_chain
      use mparameters_monomer
      use mKaps

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
      integer midR, iZ, iR, ii, zmax, zmin
       real*8 zpos, rpos, slope, ord, rcurve
      integer rpos_int

      real*8 a_param, b_param, c_param

        logical   outside
        external  outside

      real*8 x(3)

      ii = rank+1

      select case (poretype)

!-------------------- Simplest case: long pore ---------------------
      case (1)

       zpos = zposition(ii)
       rpos = (dfloat(CdimR)-0.5)*delta

       posgraft(1) = rpos
       posgraft(2) = zpos

!-------------------- empty ---------------------
      case (0)

       posgraft(1) = dfloat(dimR)/2*delta
       posgraft(2) = dfloat(dimZ)/2*delta

!-------------------- conic pore  ---------------------

      case (4)

       slope =
     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
       ord = TdimR*delta-RdimZ*delta*slope
c
c       print*, ord, slope

       zpos = zposition(ii) 
       rcurve = slope*zpos + ord ! position on the curve

       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice

       
       do while (rcurve.lt.rpos)
       rpos = rpos - delta ! is lattice site part of the membrane? 
       enddo 

       rpos = rpos + 1e-5 ! border of the lattice site

       posgraft(1) = rpos 
       posgraft(2) = zpos

!-------------------- hour glass half circle  ---------------------
       case (10)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
!       rpos = dimR*delta - delta/2.0 
       rpos = dimR*delta - 2.0*delta - 2.0*float(Kapd)*delta

       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos

       do while (outside(x))
       rpos = rpos - delta/2.0
       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos
       enddo

       posgraft(1) = rpos
       posgraft(2) = zpos

!------------------------ Hour glass like  ----------------------

       case (11)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
       rpos = dimR*delta - delta/2.0 - 2*Kapd*delta
       
       x(1) = rpos 
       x(2) = 0.0
       x(3) = zpos

       do while (outside(x))
       rpos = rpos - delta/2.0
       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos
       enddo 

       posgraft(1) = rpos
       posgraft(2) = zpos

!------------------------ cylinder  ----------------------

       case (12)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
       rpos = (dfloat(CdimR)-0.5)*delta

       posgraft(1) = rpos
       posgraft(2) = zpos

!-------------------- hour glass parabola an protein  ---------------------
c       case (5)
c
c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       zpos = zposition(ii) 
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
c
c       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice
c
c       if(rcurve.lt.rpos)rpos = rpos - delta ! is lattice site part of the membrane?
c
c       rpos = rpos - delta/2 + 1e-5 ! border of the lattice site
c
c       posgraft(1) = rpos
c       posgraft(2) = zpos


!-------------------- hour glass parabola  ---------------------
c      case (3)

c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       zpos = zposition(ii) 
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
c
c       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice
c
c       if(rcurve.lt.rpos)rpos = rpos - delta ! is lattice site part of the membrane? 
c
c       rpos = rpos - delta/2 + 1e-5 ! border of the lattice site
c
cc       posgraft(1) = rpos
c       posgraft(2) = zpos
c
       endselect
c
      if (N_chains.ne.size) then
      print*, 'The number of chains,', N_chains, ', and the
     & the number of processors',size,
     & 'should be the same'
      call MPI_FINALIZE(ierr) ! finaliza MPI
       stop
       endif

      print*, 'Chain #', rank+1, ' grafted at x=', posgraft(1), ', y=',
     & posgraft(2)

      end              

      subroutine lookup_kap

      use mlookup
      use mparameters
      use mKaps
      use mncells
      use mporesystem
      use mparameters_chain
      use mparameters_monomer

      implicit none
      integer iC, iR, iZ,  iiR, iiZ
      integer flagwall
      integer iiiZ

! define Kap

      nKap = 0
      Kaplist_cell = 0
      Kaplist_value = 0
      kapmask = 1


      do iC = 1, ncells ! loop over cells

      iR = indexa(iC,1) ! position of cell in the system
      iZ = indexa(iC,2)

! integration

      flagwall = 0  

      do iiR = 1, dimR+1 ! loop over R-neighbors
      do iiZ = -(Kapd-1)/2,(Kapd-1)/2 ! loop over Z-neighbors

      if(KapsV(iR,iiR,iiZ).ne.0d0) then ! part of Kap

      iiiZ = iiZ+iZ 

      do while (iiiZ.gt.dimZ)
      iiiZ = iiiZ-dimZ
      enddo
      do while (iiiZ.lt.1)
      iiiZ = iiiZ+dimZ
      enddo

!      if ((iiiZ.le.0).or.(iiiZ.gt.dimZ).or.(iiR.eq.dimR+1)) then ! bulk
      if (iiR.eq.dimR+1) then ! bulk

        nKap(iC) = nKap(iC) + 1
        if(nKap(iC).gt.size(Kaplist_cell,2)) then
        print*, size(Kaplist_cell,2)
        print*, 'Error in lookup_kap, increase Kaplist_cell dimensions'
        stop
        endif
        Kaplist_cell(iC,nKap(iC)) = ncells+1
        Kaplist_value(iC,nKap(iC),:) = Kaps(iR,iiR,iiZ,:)
      
      else ! not bulk
      
      if(matriz(iiR, iiiZ).gt.0)
     & then  ! inside the system and different from zero
        nKap(iC) = nKap(iC) + 1
        if(nKap(iC).gt.size(Kaplist_cell,2)) then
        print*, size(Kaplist_cell,2)
        print*, 'Error in lookup_kap, increase Kaplist_cell dimensions'
        stop
        endif
        Kaplist_cell(iC,nKap(iC)) = matriz(iiR,iiiZ)
        Kaplist_value(iC,nKap(iC),:) = Kaps(iR,iiR,iiZ,:)
      endif
      
      if(matriz(iiR,iiiZ).eq.-2) then ! part of the Kap hits the walls
      flagwall = 1
      endif

      endif ! bulk
      endif ! part of Kap

      enddo ! iiZ
      enddo ! iiR

      if(flagwall.eq.1)nKap(iC)=0 ! the Kap hit the wall 


! define mask

       if(nKap(iC).ne.0) then 
       if ((dimR-iR).le.((Kapd-1)/2))  then
       kapmask(iC) = 0
       endif
       endif

      enddo ! iC
      print*, 'Max number of neighbors', maxval(nKap)
      return
      end

      subroutine lookup_kap_s

      use mlookup
      use mparameters
      use mKaps_s
      use mKaps
      use mncells
      use mporesystem
      use mparameters_chain
      use mparameters_monomer

      implicit none
      integer iC, iR, iZ,  iiR, iiZ
      integer flagwall
      integer iiiZ

      real*8 vectr, vectz

! define Kap

      nKap_s = 0
      Kaplist_s_cell = 0
      Kaplist_s_value = 0

      do iC = 1, ncells ! loop over cells

      iR = indexa(iC,1) ! position of cell in the system
      iZ = indexa(iC,2)

! integration

      flagwall = 0  

      do iiR = 1, dimR+1 ! loop over R-neighbors
      do iiZ = -(Kapd-1)/2,(Kapd-1)/2 ! loop over Z-neighbors

      if(KapsV_s(iR,iiR,iiZ).ne.0d0) then ! part of Kap

       vectr = dfloat(iiR-iR)*delta
       vectz = dfloat(iiZ)*delta

      iiiZ = iiZ+iZ 


      do while (iiiZ.gt.dimZ)
      iiiZ = iiiZ-dimZ
      enddo
      do while (iiiZ.lt.1)
      iiiZ = iiiZ+dimZ
      enddo


!      if ((iiiZ.le.0).or.(iiiZ.gt.dimZ).or.(iiR.eq.dimR+1)) then ! bulk
      if (iiR.eq.dimR+1) then ! bulk

        nKap_s(iC) = nKap_s(iC) + 1
        if(nKap_s(iC).gt.size(Kaplist_s_cell,2)) then
        print*, size(Kaplist_s_cell,2)
        print*,'Error in lookup_kap, increase Kaplist_s_cell dimensions'
        stop
        endif
        Kaplist_s_cell(iC,nKap_s(iC)) = ncells+1
        Kaplist_s_value(iC,nKap_s(iC)) = KapsV_s(iR,iiR,iiZ)
        v0r(iC,nKap_s(iC))=vectr
        v0z(iC,nKap_s(iC))=vectz
      
      else ! not bulk
      
      if(matriz(iiR, iiiZ).gt.0)
     & then  ! inside the system and different from zero
        nKap_s(iC) = nKap_s(iC) + 1
        if(nKap_s(iC).gt.size(Kaplist_s_cell,2)) then
        print*, size(Kaplist_s_cell,2)
        print*, 'Error in lookup_kap, increase Kaplist_cell dimensions'
        stop
        endif
        Kaplist_s_cell(iC,nKap_s(iC)) = matriz(iiR,iiiZ)
        Kaplist_s_value(iC,nKap_s(iC)) = KapsV_s(iR,iiR,iiZ)
        v0r(iC,nKap_s(iC))=vectr
        v0z(iC,nKap_s(iC))=vectz
 
      endif
      
      if(matriz(iiR,iiiZ).eq.-2) then ! part of the Kap hits the walls
      flagwall = 1
      endif

      endif ! bulk
      endif ! part of Kap

      enddo ! iiZ
      enddo ! iiR

      if(flagwall.eq.1)nKap_s(iC)=0 ! the Kap hit the wall 

! define mask

      enddo ! iC
      return
      end



        subroutine create_protein

        use mprotein
        use mvariables
        use mncells
        use mlookup
        use mparameters
        use mncells
        use mporesystem

        implicit none

      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries


        integer iC, iR, iZ
        real*8 volprot, volprot2
        real*8 posz, posr, rcurve
        integer temp(0:dimR+1, 0:dimZ+1)
        real*8 protcutoff
        integer flag
        real*8 totalvol

        real*8 PdimRr
        real*8 PdimHr
        real*8 Hr
        real*8 Pposr
        real*8 center
        real*8 temp1
        integer cc, ccc
        character*5  title

        Hr = dfloat(hrange)*delta
        PdimRr = dfloat(PdimR)*delta
        Pposr = dfloat(Ppos)*delta
        PdimHr = dfloat(PdimH)*delta
       
        temp = 0

        volprot = 0.0
        volprot2 = 0.0
        protein = 0.0
        proteinC = 0.0 
        proteinq = 0.0
        proteinqC = 0.0 
        proteinhC = 0.0 
        proteinh = 0.0 
 
! assigns volume

        totalvol = 0.0

       do iZ = 1, dimZ
       do iR = 1, dimR
 
       posz = (dfloat(iZ)-0.5)*delta
       posr = (dfloat(iR)-0.5)*delta

       rcurve = 0.0 

      if((Ptype.eq.1).or.(Ptype.eq.3)) then ! sphere
      if((posz.gt.(Pposr-PdimRr))
     & .and.(posz.lt.(Pposr+PdimRr))) then
       rcurve = sqrt((PdimRr)**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.4) then ! ellipsoid
      if((posz.gt.(Pposr-PdimHr))
     & .and.(posz.lt.(Pposr+PdimHr))) then
       rcurve = PdimRr/PdimHr*sqrt((PdimHr)**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.2) then ! rounded cylinder

      temp1 = PdimRr+PdimHr/2

       if((posz.gt.(Pposr-temp1)) ! lower cap
     & .and.(posz.lt.(Pposr-PdimHr/2.0))) then
       center =  Pposr-PdimHr/2.0
       rcurve = sqrt((PdimRr)**2
     &- (posz-center)**2)
       endif


       if((posz.lt.(Pposr+temp1)) ! upper cap
     & .and.(posz.gt.(Pposr+PdimHr/2.0))) then
       center = Pposr+PdimHr/2.0
       rcurve = sqrt((PdimRr)**2
     &- (posz-center)**2)
       endif

       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
       rcurve = PdimRr 
       endif
       endif

      if(Ptype.eq.5) then ! cylinder
      temp1 = PdimRr+PdimHr/2
       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
       rcurve = PdimRr 
       endif
       endif

 

      if(rcurve.gt.posr) then ! is protein
      if(matriz(iR, iZ).eq.-3) then
      print*, 'There is a collision 
     & between the protein and the pore walls!'
      endif

      protein(iR,iZ) = prot_vol
      temp(iR,iZ) = 1
      totalvol = totalvol 
     & + (dfloat(iR)-0.5)*delta*2*pi*delta**2
      endif

      enddo
      enddo

      if(rank.eq.0)print*, 'Total prot vol', totalvol

! assigns charge to the surface

      if((Ptype.eq.1).or.(Ptype.eq.2).or.
     & (Ptype.eq.4).or.(Ptype.eq.5)) then
      do iR = 1, dimR
      do iZ = 1, dimZ
 
      flag = 0

      if(temp(iR,iZ).eq.1) then

      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1

      if(flag.eq.1) then
      volprot = volprot + delta**3*2.0*pi*(dfloat(iR)-0.5)
      proteinq(iR,iZ) = prot_q
      endif
      endif

      enddo
      enddo

      if(volprot.ne.0.0)proteinq = proteinq / volprot
      endif


      if(Ptype.eq.3) then
      do iR = 1, dimR
      do iZ = 1, dimZ

      flag = 0
      if(temp(iR,iZ).eq.1) then
      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if(flag.eq.1) then
      posz = (dfloat(iZ)-0.5)*delta
      if(posz.gt.Pposr) then
      volprot = volprot + delta**3*2.0*pi*(dfloat(iR)-0.5)
      else
      volprot2 = volprot2 + delta**3*2.0*pi*(dfloat(iR)-0.5)
      endif
      endif
      endif
      enddo
      enddo


      do iR = 1, dimR
      do iZ = 1, dimZ
      flag = 0
      if(temp(iR,iZ).eq.1) then
      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if(flag.eq.1) then
       posz = (dfloat(iZ)-0.5)*delta
      if(posz.gt.Pposr) then
      proteinq(iR,iZ) = prot_q/volprot
      else
      proteinq(iR,iZ) = prot_q2/volprot2
      endif
      endif
      endif
      enddo
      enddo

      end if



! assigns hydrophobic interaction matrix

        do iZ = 1, dimZ
        do iR = 1, dimR

       posz = (dfloat(iZ)-0.5)*delta
       posr = (dfloat(iR)-0.5)*delta

       rcurve = 0.0

      if((Ptype.eq.1).or.(Ptype.eq.3)) then ! sphere
       if((posz.gt.(Pposr-PdimRr-Hr))
     & .and.(posz.lt.(Pposr+PdimRr+Hr))) then
       rcurve = sqrt(((PdimRr+Hr))**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.2) then ! rounded cylinder

      temp1 = PdimRr+PdimHr/2

       if((posz.gt.(Pposr-temp1-Hr)) ! lower cap
     & .and.(posz.lt.(Pposr-PdimHr/2.0))) then
       center =  Pposr-PdimHr/2.0
       rcurve = sqrt((PdimRr+Hr)**2
     &- (posz-center)**2)
       endif


       if((posz.lt.(Pposr+temp1+Hr)) ! upper cap
     & .and.(posz.gt.(Pposr+PdimHr/2.0))) then
       center = Pposr+PdimHr/2.0
       rcurve = sqrt((PdimRr+Hr)**2
     &- (posz-center)**2)
       endif

       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
      rcurve = PdimRr+Hr
       endif
       endif

      if(Ptype.eq.5) then ! cylinder
      temp1 = PdimRr+PdimHr/2
       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
      rcurve = PdimRr+Hr
       endif
       endif


      if(Ptype.eq.4) then ! ellipsoid
      if((posz.gt.(Pposr-PdimHr))
     & .and.(posz.lt.(Pposr+PdimHr))) then
       rcurve = PdimRr/PdimHr*sqrt((PdimHr)**2
     &- (posz-Pposr)**2)
      endif
      endif


      if(rcurve.gt.posr) then ! is protein
      if((matriz(iR, iZ).ge.1).and.(protein(iR, iZ).eq.0.0)) then ! if the point is inside the system and is not protein 
      proteinh(iR, iZ) = 1.0
      endif
      endif

      enddo
      enddo

! create proteinC and proteinqC

      do iC = 1, ncells
      proteinC(iC) = protein(indexa(iC,1),indexa(iC,2))
      proteinqC(iC) = proteinq(indexa(iC,1),indexa(iC,2))
      proteinhC(iC) = proteinh(indexa(iC,1),indexa(iC,2))
      enddo

      if(savetodisk_flag.eq.12) then

      if(rank.eq.0) then
      print*, 'Saving protein properties and exiting'
      cc = 1
      ccc = 1

      title = 'protC'
      call savetodisk(proteinC, title, cc ,ccc)

      title = 'proqC'
      call savetodisk(proteinqC, title, cc ,ccc)

      title = 'prohC'
      call savetodisk(proteinhC, title, cc ,ccc)
      endif
       
      call MPI_FINALIZE(ierr) ! finalize MPI
      stop
      endif 

      end

