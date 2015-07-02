      subroutine fkfun(x,f,ier)

      use mlookup
      use mchains
      use mparameters
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mparameters_monomer
      use mprotein
      use mKaps
      use mKaps_s

!
! fbound(im,iC,ii) = fraction of segments of type im and located at pC that are bound to a kap with center in iC
!                    pC is the position of the ii receptor of the kap located at iC
! funbound(im, pC) = fraction of the segments of type im and located at pC that are free
! recetor(iC) = total concentration of receptors at iC
!
!
!



      implicit none
 
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      double precision Factorcurv

      integer ini, fin

      integer ntot
      real*8 x((2+N_poorsol)*ncells)
      real*8 f((2+N_poorsol)*ncells)

      real*8  protemp, protemp1

      integer i,j, k, ix, iy, iC, ii, aR, temp, aZ, jj, im, imm

      real*8 temp2

      real*8 psi(ncells+2) ! psi se define asi para evitar problemas al construir las fs
      real*8 xtotal(N_poorsol, ncells+2) ! psi se define asi para evitar problemas al construir las fs

      real*8 xh(ncells) ! psi se define asi para evitar problemas al construir las fs
 
! Kinsol

      integer*8 neq
      integer*4 ier

! MPI

       integer spp ! sitios por procesador

       integer cuantas_p(2)

! solving auxiliary variables

       real*8 avpol_temp(N_monomer, N_chains+2, ncells+2)
       real*8 avpolsum(N_monomer, dimR*dimZ)
       real*8 receptor_temp(ncells+2)
       real*8 avpol_tosend(N_monomer, N_chains+2, ncells)

       real*8 xpot(N_monomer,ncells+2)
       real*8 xpotp(N_monomer,ncells+2)
     
! bit operations

       integer*1 displ_one(0:7)
       integer*4 pC, pZ, pR, pR2

      integer*1 displacement(maxlong+1,2)
      integer*1 binary(int(maxlong/2))

! kap 
       real*8 prokap
       real*8 tempprot

       integer pC2

C-----------------------------------------------------

! Jefe

      if(rank.eq.0) then ! llama a subordinados y pasa vector x
       flagsolver = 1

       CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)

       CALL MPI_BCAST(x, (2+N_poorsol)*ncells , MPI_DOUBLE_PRECISION,0
     &  , MPI_COMM_WORLD,err)

      endif


! Recupera xh y psi

      ntot = ncells

            do iC = 1,ntot

            xh(iC)=x(iC)
            psi(iC)=x(iC+ntot)

            do i=1, N_poorsol
            xtotal(i,iC)=x(iC+(1+i)*ntot) ! v-frac of the good solvent polymers
            enddo

            enddo

! Condiciones de borde psi 
 
! This implictly considers SIGMAQ = 0
! 


! bulk
      psi(ncells+1) = 0.0
      psi(ncells+2) = 0.0 

! lookuptable has teh boundary conditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Volume fraction of mobile species and dissociation fraction

           avpol=0.0
           avpol_tosend = 0.0
           avpol_temp = 0.0
           fdis = 0.0
           fdisP = 0.0
           rhokap = 0.0
           rhokapb = 0.0
           receptor = 0.0
           receptor_temp = 0.0 
           fbound = 0.0
           funbound = 0.0

           do iC = 1, ncells
     
           xpos(iC) = expmupos*(xh(iC)**vsalt)
     &     *dexp(-psi(iC)*zpos) ! ion plus volume fraction

           xneg(iC) = expmuneg*(xh(iC)**vsalt)
     &     *dexp(-psi(iC)*zneg) ! ion neg volume fraction
     
           xHplus(iC) = expmuHplus*(xh(iC))
     &     *dexp(-psi(iC))           ! H+ volume fraction
     
           xOHmin(iC) = expmuOHmin*(xh(iC))
     &     *dexp(+psi(iC))           ! OH-  volume fraction

         if (prot_q < 0) then
               fdisP(iC) =
     &              1.0 /(1.0 + xHplus(iC)/(K0P*xh(iC)) )
         else ! base
               fdisP(iC) =
     &              1.0 /(1.0 + xOHmin(iC)/(K0P*xh(iC)) )
         endif

 
         do im =1,N_monomer
            if (zpol(im).eq.1) then !BASE
               fdis(im,iC) =
     &              1.0 /(1.0 + xOHmin(iC)/(K0(im)*xh(iC)) )
            else if (zpol(im).eq.-1) then !ACID
               fdis(im,iC) =
     &              1.0 /(1.0 + xHplus(iC)/(K0(im)*xh(iC)) )
            endif
         enddo

            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xpot boundary condition

          do i = 1, N_poorsol
          xtotal(i, ncells+1) = xtotalbulk(i)
          xtotal(i, ncells+2) = xtotalbulk(i)
          enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xpot calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do im =1, N_monomer
      do iC = 1, ntot

      xpot(im,iC) = dlog(xh(iC))*vpol

      if(zpol(im).ne.0) then
       xpot(im,iC) = xpot(im,iC)
     & -psi(iC)*zpol(im)
     & -dlog(fdis(im, iC))
      endif

      if(hydroph(im).ne.0) then 

               xpot(im,iC) = xpot(im,iC)+henergy(hydroph(im))
     & *proteinhC(iC)*hst

               

               do jj = 1, nXu(iC) ! loop over kai neighbors
               do ii = 1, N_poorsol ! loop over different poor sv types

                        xpot(im,iC) = xpot(im,iC) +
     &                       (st_matrix(hydroph(im),ii) ! st_matrix(x, y) : interaction of hydrophobic segments of type x with those of type y 
     &                       *st/(vsol*vpol)*           ! st in kT/monomer          
     &                       Xulist_value(iC,jj)*
     &                       xtotal(ii, Xulist_cell(iC, jj)))

               enddo ! ii
               enddo  !jj

      endif ! hydrophob
      enddo
      enddo

! bulk
      xpot(:,ncells+1)=xpotbulk(:)
      xpot(:,ncells+2)=xpotbulk(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation kap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       ini = int(float(rank)*float(ntot)/float(size))+1
!       fin = int(float(rank+1)*float(ntot)/float(size))
 

       if(rank.eq.0) then
       do iC = 1, ntot

       pR = indexa(iC, 1)
       pZ = indexa(iC, 2)
       
       prokap = 0.0

       tempprot = 0.0

       do jj = 1, nKap(iC) ! loop over sites with segments
         pC2 = kaplist_cell(iC,jj) ! cell number jj
         do im = 1, N_monomer ! loop over monomer types 
           prokap = prokap + xpot(im, pC2)*(kaplist_value(iC,jj,im)) ! number of aminoacids of type im that the kap in iC has in neighbor cell jj
         enddo
       enddo

       if(nKap(iC).ne.0) then

! rhokap
         rhokap(iC) = expmukap*exp(prokap)/(vkap)

! avpol
         do jj = 1, nKap(iC) ! loop over sites with aminoacids
         pR2 = indexa(kaplist_cell(iC,jj),1)
         pC2 = kaplist_cell(iC,jj)
         do im = 1, N_monomer
         avpol_temp(im,N_chains+1,pC2) =
     &   avpol_temp(im,N_chains+1,pC2)
     &   + expmukap*exp(prokap)*vpol/vkap
     &   * kaplist_value(iC,jj,im)
     &   *(dfloat(pR)-0.5)/(dfloat(pR2)-0.5)
         enddo ! im
         enddo ! jj
! receptor

         do jj = 1, nKap_s(iC)
         pR2 = indexa(kaplist_s_cell(iC,jj),1)
         pC2 = kaplist_s_cell(iC,jj)
         receptor_temp(pC2) =
     &   receptor_temp(pC2)
     &   + rhokap(iC)
     &   * kaplist_s_value(iC,jj)
     &   * (dfloat(pR)-0.5)/(dfloat(pR2)-0.5)
         enddo ! jj

       else
         rhokap(iC) = 0.0
       endif
       enddo ! loop over center



! apply mask to deal with bulk boundaries

       do iC = 1, ntot


       do im = 1, N_monomer
       avpol_temp(im,N_chains+1,iC)=
     & avpol_temp(im,N_chains+1,iC)*kapmask(iC) + 
     & xkapbulk*nmonkap(im)/sum(nmonkap)*(1.0-kapmask(iC))
       enddo

       receptor_temp(iC)=
     & receptor_temp(iC)*kapmask(iC) +
     & xkapbulk/vkap*(1.0-kapmask(iC))

       enddo
       endif

! receptor
       receptor(1:ncells) = receptor_temp(1:ncells)

! MPI: send receptor to all processors

      call MPI_BCAST(receptor
     & , dimR*dimZ,
     &   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

      CALL MPI_BARRIER(MPI_COMM_WORLD,err)

! Calculation of funbound

         do iC = 1, ncells
         do im =1,N_monomer
            if (rbind(im).eq.1) then ! RECEPTOR
               funbound(im,iC) =
     &              1.0 /(1.0 + Kbind0(im)*receptor(iC))
            endif
         enddo
         enddo

! xpot for polymer

      xpotp = xpot

      do im =1, N_monomer
      do iC = 1, ntot

      if(rbind(im).eq.1) then
       xpotp(im,iC) = xpotp(im,iC)
     & -dlog(funbound(im, iC))
      endif

      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation chains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ii = rank+1


      q = 0.0
      sumprolnpro = 0.0
      endtoend_av = 0.0

      do i=1,newcuantas(ii)

         lnpro=0.0

!! DECODER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         pC = fs(i)                                        !
!         pR = indexa(pC,1)                                 !
!         pZ = indexa(pC,2)                                 !
!         binary(:) = displ(i,:)                            !
!         call decode(displacement, binary, long(ii))                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,long(ii)

!         pR = pR + displacement(j,1) 
!         pZ = pZ + displacement(j,2) 
!         if(pZ.eq.0)pZ=dimZ
!         if(pZ.eq.(dimZ+1))pZ=1

         pC = inc(i, j)
         pR = indexa(pC,1)
         pZ = indexa(pC,2)

!         print*, i, pR, pZ

         lnpro = lnpro + xpotp(segtype(ii, j),pC)

         enddo ! j
            
            pro = dexp(lnpro)

            q=q+pro
            sumprolnpro = sumprolnpro + pro*lnpro

            endtoend_av = endtoend_av + pro*endtoend(i)

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         pC = fs(i)                                     !
!         pR = indexa(pC,1)                              !
!         pZ = indexa(pC,2)                              !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,long(ii)

!         pR = pR + displacement(j,1)
!         pZ = pZ + displacement(j,2)
!         if(pZ.eq.(dimZ+1))pZ=1
!         if(pZ.eq.0)pZ=dimZ
       
!         pC = matriz(pR, pZ)

         pC = inc(i, j)
         pR = indexa(pC,1)
         pZ = indexa(pC,2)

         avpol_temp(segtype(ii,j), ii, pC) 
     & = avpol_temp(segtype(ii,j), ii, pC)
     & + pro*chainsperdelta(ii)*vsol
     & *vpol*Factorcurv(pR)

         enddo ! j
         enddo !i

          avpol_temp(:,ii,:)=avpol_temp(:,ii,:)/q
          endtoend_av = endtoend_av/q

          avpol_tosend(:,:,1:ncells) = avpol_temp(:,:,1:ncells)

c---------------------------------------------------------------------        
c------------------ MPI ----------------------------------------------
c---------------------------------------------------------------------

         call MPI_REDUCE(avpol_tosend, avpol
     & , ncells*N_monomer*(N_chains+2),
     &   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

      CALL MPI_BARRIER(MPI_COMM_WORLD,err)

      if(rank.ne.0)goto 3333
!!!!!!!!!!! IMPORTANT, SLAVES STOP HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c ------------------------------------------------------------------------
c Calculate volume and charge of bound proteins 
c ------------------------------------------------------------------------

      avpolsum = 0.0
      do ii = 1, N_chains
      avpolsum(:,:) = avpolsum(:,:) + avpol(:,ii,:)
      enddo

      do imm = 1, N_monomer

      if(rbind(imm).eq.1) then ! this monomer type has binding

      do iC = 1, ncells ! center of Kap

      if(nKap(iC).ne.0) then

         do ii = 1, nKap_s(iC) ! loop over the sites of the Kap that have receptors
 
         pC = kaplist_s_cell(iC,ii) ! cell of the receptor
         pR = indexa(pC, 1) ! r coordinate of the receptor

         pR2 = indexa(iC, 1) ! r coordinate of the particle

         fbound(imm,iC,ii)=
     & funbound(imm,pC)*Kbind0(imm)*kaplist_s_value(iC,ii)*rhokap(iC)
     & *(dfloat(pR2)-0.5)/(dfloat(pR)-0.5)

! add avpol

         do jj = 1, nKap(iC) ! loop over sites with aminoacids

         pR2 = indexa(kaplist_cell(iC,jj),1)
         pC2 = kaplist_cell(iC,jj)

         do im = 1, N_monomer
         avpol(im,N_chains+2,pC2) =
     &   avpol(im,N_chains+2,pC2)
     &   + avpolsum(imm,pC)/vpol
     &   *fbound(imm,iC,ii)*kaplist_value(iC,jj,im)*vpol
     &   *(dfloat(pR)-0.5)/(dfloat(pR2)-0.5)
         enddo ! im
         enddo ! jj

         pR2 = indexa(iC,1)

          rhokapb(iC) = rhokapb(iC) + avpolsum(imm,pC)/vpol
     &   *fbound(imm,iC,ii)
     &   *(dfloat(pR)-0.5)/(dfloat(pR2)-0.5)

         enddo ! Ii
      endif
      enddo ! iC 

      endif
      enddo ! imm

C----------------------------------------------------------------------------------------------
C   Construct f vector
C----------------------------------------------------------------------------------------------

! Qtot

      do iC=1,ntot

         qtot(iC) = 
     &   (zpos*xpos(iC)+zneg*xneg(iC))/vsalt 
     &   + xHplus(iC)-xOHmin(iC)

      if(weakP.eq.1) then
      qtot(iC) = qtot(iC)  + proteinqC(iC)*vsol*fdisP(iC) 
      else
      qtot(iC) = qtot(iC)  + proteinqC(iC)*vsol
      endif

      do im = 1, N_monomer
      do ii = 1, N_chains+2

      qtot(iC) = qtot(iC) + avpol(im, ii, iC)*fdis(im, iC)*zpol(im)/vpol

      enddo ! ii
      enddo ! im
 
      enddo

! Volume fraction

      do iC=1, ntot
         
               f(iC)=
     &   xh(iC) + xneg(iC) 
     &  + xpos(iC) + xHplus(iC) + xOHmin(iC)
     &  + proteinC(iC) - 1.0d0
   
      do im = 1, N_monomer
      do ii = 1, N_chains+2

      f(iC) = f(iC) + avpol(im, ii, iC)

      enddo ! ii
      enddo ! im

      enddo

! Poisson eq.

            do iC=1,ntot

! Cilindro (derivada al centro), ver notas     
     
               f(iC+ntot)=
     & psi(rp(iC)) -2*psi(iC) + psi(rm(iC)) +
     & (0.5/(dfloat(indexa(iC,1))-0.5))*(psi(rp(iC))-psi(rm(iC))) + ! termino de curvatura
     & psi(zp(iC)) -2*psi(iC) + psi(zm(iC)) + ! derivada en z
     & qtot(iC)*constq
     
               f(iC+ntot)=f(iC+ntot)/(-2.0) ! mejora kinsol...
      enddo

! poor solvent

      do ii = 1, N_poorsol
      do iC = 1, ntot

      f(iC+(1+ii)*ntot) = xtotal(ii,iC)

      do im = 1, N_monomer
      if(hydroph(im).eq.ii) then
      do jj = 1, N_chains+2
        f(iC+(1+ii)*ntot) = f(iC+(1+ii)*ntot) - avpol(im,jj,iC)
      enddo ! jj
      endif
      enddo ! im
      enddo ! iC
      enddo ! ii

      iter = iter + 1

      norma = 0.0

      do i = 1, (2+N_poorsol)*ntot
      norma = norma +(f(i))**2    
      enddo
      
      print*, iter, norma, q
      open(unit=6660, file='iter.dat', access='append')
      write(6660,*)iter, norma, q
      close (6660)
      ier = 0.0

! saves infile

      if(mod(iter, save_every).eq.0) then
      print*, 'Save resume file'
      open(unit=45, file = 'out.temp.dat')
      do i = 1, (2+N_poorsol)*ncells
      write(45, *)x(i)
      enddo
      close(45)
      endif



 3333 ier = 0.0

      return
      end


