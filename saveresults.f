      subroutine saveresults(cc, ccc) ! save results

      use mlookup
      use mparameters ! parametros     
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mprotein
      use mparameters_monomer
      use mKaps_s
      use mKaps

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries


      real*8 movpos,movneg,movHplus,movOHmin,Gporo,Gvacio
      real*8 radio
      integer iZ, pZ
      real*8 avvect(dimR*dimZ,2)
      real*8 oop(dimR*dimZ)
      real*8 oops(dimR*dimZ)
      real*8 avpol_monom_pol(N_monomer, dimR*dimZ)
      integer ccc, i, ii, jj, j, cc, im, is,iC,iR
      real*8 avpol_all(dimR*dimZ)
      real*8 avpol_temp(dimR*dimZ)
      real*8 fdis_temp(dimR*dimZ)
      character*5  title
      character*24 filename 
      real*8 temp
      real*8 temp1

      real*8 sumkaprho
      real*8 avpolsum(N_monomer, dimR*dimZ)
      integer pC, pR

C----------------------------------------------------------
C  OUTPUT!
C----------------------------------------------------------

      print*, 'SAVETODISK', savetodisk_flag

      if(rank.eq.0) then ! only master can save to disk

!!!!!!!!!!!!!!!!!!! saves files  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! POLYMER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Polymer, sum

      if(savetodisk_flag.ne.10) then
      avpol_all(:) = 0
      do im = 1, N_monomer
      do ii = 1, N_chains
      avpol_all(:) = avpol_all(:) + avpol(im, ii, :)
      enddo
      enddo
      title = 'avpol'
      call savetodisk(avpol_all, title, cc ,ccc)
      endif


! Fraction hydrophobic - hydrophylic
      if((savetodisk_flag.eq.3).or.savetodisk_flag.eq.0) then
      do j = 1, N_poorsol
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(hydroph(im).eq.j) then 
         do ii = 1, N_chains+2
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo
      endif
      enddo ! im
      write(title,'(A4, I1.1)')'hpho',j
      call savetodisk(avpol_temp, title, cc ,ccc)
      enddo ! j

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(hydroph(im).eq.0) then
         do ii = 1, N_chains+2
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo
      endif 
      enddo ! im
      write(title,'(A5)')'hphil'
      call savetodisk(avpol_temp, title, cc ,ccc)

! pos, neg, neutral
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.0) then
         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo
      endif
      enddo ! im
      write(title,'(A5)')'neutr'
      call savetodisk(avpol_temp, title, cc ,ccc)

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.1) then
         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo
      endif
      enddo ! im
      write(title,'(A5)')'posit'
      call savetodisk(avpol_temp, title, cc ,ccc)

!
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.-1) then
         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo
      endif
      enddo ! im
      write(title,'(A5)')'negat'
      call savetodisk(avpol_temp, title, cc ,ccc)

      endif

! all polymers
      if((savetodisk_flag.eq.2).or.(savetodisk_flag.eq.0)) then
      do ii = 1, N_chains
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      enddo
      write(title,'(A3, I2.2)')'cha',ii
      call savetodisk(avpol_temp, title, cc ,ccc)
      enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! kap, sum

      if(savetodisk_flag.ne.10) then


      avpol_all(:) = 0
      do im = 1, N_monomer
      ii = N_chains+1
      avpol_all(:) = avpol_all(:) + avpol(im, ii, :)
      enddo
      title = 'avkap'
      call savetodisk(avpol_all, title, cc ,ccc)

      avpol_all(:) = 0
      do im = 1, N_monomer
      ii = N_chains+2
      avpol_all(:) = avpol_all(:) + avpol(im, ii, :)
      enddo
      title = 'avkbn'
      call savetodisk(avpol_all, title, cc ,ccc)
      endif



! Kap number density

      if(savetodisk_flag.ne.10) then
      title = 'rhokp'
      call savetodisk(rhokap, title, cc ,ccc)

      title = 'rhokb'
      call savetodisk(rhokapb, title, cc ,ccc)
      endif



      if(savetodisk_flag.ne.10) then
      title = 'recep'
      call savetodisk(receptor, title, cc ,ccc)
      endif


! vectors

! avpol preliminary calculations

      if(savetodisk_flag.ne.10) then

      avpol_monom_pol = 0.0
      do im = 1, N_monomer
      do ii = 1, N_chains
      avpol_monom_pol(im, :) = avpol_monom_pol(im, :) +  avpol(im,ii,:)
      enddo
      enddo

      do im = 1, N_monomer

       

      avvect = 0.0
      oop = 0.0
      oops = 0.0
      do iC = 1, ncells
      do ii = 1, nKap_s(iC)

      pC = kaplist_s_cell(iC,ii)

      if(avpol_monom_pol(im, pC).ne.0.0) then

      iR = indexa(iC,1)
      iZ = indexa(iC,2)
      pR = indexa(pC,1)
      pZ = indexa(pC,2)

      avvect(pC,1) = avvect(pC,1) + 
     & avpol_monom_pol(im, pC)*fbound(im,iC,ii)
     & *v0r(iC,ii)

      avvect(pC,2) = avvect(pC,2) + 
     & avpol_monom_pol(im, pC)*fbound(im,iC,ii)
     & *v0z(iC,ii)

      if((v0r(iC,ii)**2 + v0z(iC,ii)**2).ne.0.0) 
     & then
      oop(pC) = oop(pC) + ((3.0*(v0r(iC,ii))**2
     & /((v0r(iC,ii))**2 + (v0z(iC,ii)**2)))-1.0)/2.0*
     & fbound(im,iC,ii)
      endif

      oops(pC) = oops(pC)+fbound(im,iC,ii)

      endif

      enddo ! ii
      enddo ! iC
   
      do iC = 1, ncells
      if(oops(iC).ne.0.0)oop(iC)=oop(iC)/oops(iC)
      enddo

      write(title,'(A3, I2.2)')'vec',im
      call savevect(avvect, title, cc ,ccc)

      write(title,'(A3, I2.2)')'oop',im
      call savetodisk(oop, title, cc ,ccc)
 
      enddo ! im
 
      endif



! Fraction hydrophobic - hydrophylic

      if((savetodisk_flag.eq.3).or.savetodisk_flag.eq.0) then

      do j = 1, N_poorsol
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(hydroph(im).eq.j) then 
         ii = 1+N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      endif
      enddo ! im
      write(title,'(A4, I1.1)')'phok',j
      call savetodisk(avpol_temp, title, cc ,ccc)
      enddo ! j
       
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(hydroph(im).eq.0) then
         ii =  1+N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      endif 
      enddo ! im
      write(title,'(A5)')'philk'
      call savetodisk(avpol_temp, title, cc ,ccc)

! pos, neg, neutral

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.0) then
         ii = 1+N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      endif
      enddo ! im
      write(title,'(A5)')'neutk'
      call savetodisk(avpol_temp, title, cc ,ccc)
!
      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.1) then
         ii = 1+N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      endif
      enddo ! im
      write(title,'(A5)')'positk'
      call savetodisk(avpol_temp, title, cc ,ccc)

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.-1) then
         ii = 1+N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
      endif
      enddo ! im
      write(title,'(A5)')'negatk'
      call savetodisk(avpol_temp, title, cc ,ccc)

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Other species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if((savetodisk_flag.eq.0).or.(savetodisk_flag.eq.5)) then

! Total charge
      title = 'qtodo'
      call savetodisk(qtot, title, cc ,ccc)

! Solvente
      title = 'avsol'
      call savetodisk(xsol, title, cc, ccc)

! Cationes
      title = 'avpos'
      call savetodisk(xpos, title, cc, ccc)

! Aniones
      title = 'avneg'
      call savetodisk(xneg, title, cc, ccc)

! H+
      title = 'avHpl'
      call savetodisk(xHplus, title, cc, ccc)

! OH-
      title = 'avOHm'
      call savetodisk(xOHmin, title, cc,ccc)

! fdis
      fdis_temp(:) = 0

      do im = 1, N_monomer

      fdis_temp(:) = fdis(im,:)

      write(title,'(A3, I1.1)')'fdis',im
      call savetodisk(fdis_temp, title, cc ,ccc)

      enddo


! funbound

      fdis_temp(:) = 0
      do im = 1, N_monomer
      fdis_temp(:) = funbound(im,:)
      write(title,'(A3, I1.1)')'fun',im
      call savetodisk(fdis_temp, title, cc ,ccc)
      enddo

! funbound_diff

      do im = 1, N_monomer
      fdis_temp(:) = 1.0
      do iC = 1, ncells

      do ii = 1, nKap_s(iC)
      pC = kaplist_s_cell(iC,ii)
      fdis_temp(pC) = fdis_temp(pC)-fbound(im,iC,ii)
      enddo ! ii

      enddo ! iC

      write(title,'(A3, I1.1)')'fdf',im
      call savetodisk(fdis_temp, title, cc ,ccc)
      enddo

! fdisP
      title = 'fdisP_'
      call savetodisk(fdisP, title, cc ,ccc)

      endif ! savetodisk flag

       if(((savetodisk_flag.eq.0).or.(savetodisk_flag.eq.5)
     & .or.(savetodisk_flag.eq.11))) then
! protein
      title = 'protC'
      call savetodisk(proteinC, title, cc ,ccc)

! Potencial electrostatico
      title = 'poten'
      call savetodisk(psi2, title, cc, ccc)
      endif

!!!!!!!!!!!!!!!!! Informacion del sistema !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         write(filename,'(A7, I3.3, A1, I3.3, A4)')'system.',ccc,'.',cc,
     &     '.dat'
         open (unit=510, file=filename)

         write(510,*)'st          = ',st ! residual size of iteration vector
         write(510,*)'fnorm       = ',norma ! residual size of iteration vector
         write(510,*)'delta       = ',delta
         write(510,*)'vsol        = ',vsol
         write(510,*)'vsol        = ',vpol
         write(510,*)'vsalt       = ',vsalt*vsol
         write(510,*)'csalt       = ',csalt
         write(510,*)'pHbulk      = ',pHbulk
         write(510,*)'pKw         = ',pKw
         write(510,*)'zpos        = ',zpos
         write(510,*)'zneg        = ',zneg
         write(510,*)'cuantas     = ',cuantas
         write(510,*)'newcuantas     = ',newcuantas
         write(510,*)'iterations  = ',iter
         write(510,*)'mukap       = ',dlog(expmukap)
         write(510,*)'CHAINS PER DELTA  = ',chainsperdelta
         write(510,*)'system length / nm = ', dimZ*delta
         write(510,*)'sites with grafted chains = ', CdimZ
         write(510,*)'Kbind0 = ', Kbind0
         write(510,*)'mukap = ', log(expmukap)

!!!!!!!!!!!!!!!!!!!!!! CHECK TOTAL NUMBER OF CHAINS !!!!!!!!!!!!!!!!!
        temp = 0.0
        
        do iC = 1, ncells
        do im = 1, N_monomer
        do ii = 1, N_chains

        iR = indexa(iC, 1)
        temp = temp + avpol(im,ii,iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2

        enddo
        enddo
        enddo

        temp = temp/vpol/vsol
        write(510,*)'Total number of segments in system', temp

        temp = 0.0
        do ii = 1, N_chains
        temp = temp + long(ii)*chainsperdelta(ii)
        enddo
        write(510,*)'Total number of segments in system should be',
     &  temp

!!!!!!!!!!!!!!!!!!!!!!!! CHECK TOTAL NUMBER OF BOUND PROTEINS !!!!!!!!!

! 1: from avpol


        avpolsum = 0.0
        do ii = 1, N_chains
        avpolsum(:,:) = avpolsum(:,:) + avpol(:,ii,:)
        enddo

        temp = 0.0

        do iC = 1, ncells ! center of the kap
 
        iR = indexa(iC,1)

        do im = 1, N_monomer
        if(rbind(im).eq.1) then

        do ii = 1, nKap_s(iC)

        pC = kaplist_s_cell(iC,ii)
        pR = indexa(pC,1)

        temp = temp + fbound(im,iC,ii)*avpolsum(im,pC)/vpol*
     &  (dfloat(pR)-0.5)/(dfloat(iR)-0.5)*
     &  (dfloat(iR)-0.5)*delta*2.0*pi*delta**2

        enddo ! ii
        endif ! im
        enddo ! im

        enddo ! iC

        write(510,*)'Total number of adsorbed proteins (from avpol)',
     &  temp/vsol


! 2: from rhokapb

        temp = 0.0

        do iC = 1, ncells
        iR = indexa(iC,1)
        temp = temp + rhokapb(iC)*
     &  (dfloat(iR)-0.5)*delta*2.0*pi*delta**2
        enddo

        write(510,*)'Total number of adsorbed proteins (from rhokapb)',
     &  temp/vsol

! 3: from volumen

        temp = 0.0

        do iC = 1, ncells
        iR = indexa(iC,1)
        do im = 1, N_monomer
        temp = temp + avpol(im,N_chains+2,iC)/vkap*
     &  (dfloat(iR)-0.5)*delta*2.0*pi*delta**2
        enddo
        enddo

        write(510,*)'Total number of adsorbed proteins (from vkap)',
     &  temp/vsol


! Only inside pore


! from volumen

        temp = 0.0
        temp1 = 0.0

        do iC = 1, ncells
        iR = indexa(iC,1)
        iZ = indexa(iC,2)

        if((iZ.gt.((dimZ-CdimZ)/2))
     &  .and.(iZ.le.((dimZ-CdimZ)/2+CdimZ))) then

        do im = 1, N_monomer
        temp = temp + avpol(im,N_chains+2,iC)*
     &  (dfloat(iR)-0.5)*delta*2.0*pi*delta**2

        do ii = 1, N_chains 
        temp1 = temp1 + avpol(im,ii,iC)*
     &  (dfloat(iR)-0.5)*delta*2.0*pi*delta**2
        enddo

        enddo
        endif

        enddo ! iC

        write(510,*)'%F (protein)',temp
     & /(pi*(dfloat(CdimR)**2)*CdimZ*delta**3)
        write(510,*)'%F (polymer)',temp1
     & /(pi*(dfloat(CdimR)**2)*CdimZ*delta**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PARTICLE DISOC !!!!!!!!!!!!!!!!!
        if(weakP.eq.1) then

        temp = 0.0
        temp1 = 0.0

        do iC = 1, ncells
        iR = indexa(iC, 1)
        temp = temp + proteinqC(iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2
        enddo

        write(510,*)'Total number of cargo ionizable charges', temp


        do iC = 1, ncells
         
        iR = indexa(iC, 1)
        temp1 = temp1 + proteinqC(iC)*fdisP(iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2

        enddo

        write(510,*)'Total number of cargo ionized charges', temp1

        write(510,*)'Fraction of cargo ionized charges', temp1/temp
        endif
!!!!!!!!!!!!!!!!!!!!!! NUMBER OF KAPS !!!!!!!!!!!!!!!!!
        temp = 0.0
        
        do iC = 1, ncells
        iR = indexa(iC, 1)
        temp = temp + rhokap(iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2
        enddo

        temp = temp/vpol/vsol
        write(510,*)'Total number of kaps in system', temp

!!!!!!!!!!!!!!!!!!!!!!!!!!! conductivity, only ptype = 1 !!!!!!!!


        if (poretype.eq.1) then

        radio = delta * CdimR ! radio del poro en nm
        movpos = 7.352 ! mobility S/m/M
        movneg = 7.634
        movHplus = 34.982
        movOHmin = 19.8               

       Gvacio = Gvacio + xposbulk/vsalt/vsol*movpos
       Gvacio = Gvacio + xnegbulk/vsalt/vsol*movneg
       Gvacio = Gvacio + xHplusbulk/vsol*movHplus
       Gvacio = Gvacio + xOHminbulk/vsol*movOHmin
       Gvacio = Gvacio * 1e24/Na ! Corrige unidades concentracion
       Gvacio = Gvacio * 1.0d-18/1d-6 ! Corrige unidades
       Gvacio = Gvacio * pi * radio**2 ! Corrige unidades concentracion

       if(mod(dimZ,2).eq.0) then
         iz = int(dimZ/2) 
       else
         iz = int((dimZ+1)/2) 
       endif

        Gporo = 0.0
        do iR = 1, CdimR

        iC = matriz(iR,iZ)

        Gporo=Gporo+xpos(iC)/vsalt/vsol*movpos*(dfloat(iR)-0.5)
        Gporo=Gporo+xneg(iC)/vsalt/vsol*movneg*(dfloat(iR)-0.5)
        Gporo=Gporo+xHplus(iC)/vsol*movHplus*(dfloat(iR)-0.5)
        Gporo=Gporo+xOHmin(iC)/vsol*movOHmin*(dfloat(iR)-0.5)
        enddo

        Gporo = Gporo * 1e24/Na ! Corrige unidades concentracion
        Gporo = Gporo * 1.0d-18/1d-6 ! Corrige unidades
        Gporo = Gporo * 2*pi * delta**2


        write(510,*)'Conductance for 1 um (empty)', Gvacio
        write(510,*)'Conductance for 1 um', Gporo

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         close(510)
 
         endif ! rank
        
      end
      

C**************************************************************
C*************************************************************
c
c

      subroutine savetodisk(array2, title, counter, counter2)

      use mlookup
      use mparameters
      use mncells
      use mparameters_chain
      use mparameters_monomer

       implicit none
       real*4 single
       integer ix, iy, iR, i, jx, jy, jz, iZ, iC
       integer maxT, iT
       real*8 posx, posy, posz
       character*30 filename, tempc
       character*6 titlez
       character*5 title
       integer counter, counter2
       real*8 arrayC(ncells)
       real*8 array2(dimR*dimZ)

       real*8 array(dimR, dimZ)


! Variables

! Crea array

      array = -10 ! hay que extraer las celdas con valor -1000

      do iC = 1, ncells
      array(indexa(iC,1),indexa(iC,2))=array2(iC)
      enddo

! Archivo paraview 3D

      if((savetodisk_type.eq.0).or.(savetodisk_type.eq.2)) then

      maxT = 1

      write(filename,'(A5, A1,I3.3, A1, I3.3, A4)')title,'.',
     &  counter, '.',  counter2, '.vtk'
      open (unit=45, file=filename)
      write(45,'(A)')'# vtk DataFile Version 2.0'
      write(45,'(A)')title
      write(45,'(A)')'ASCII'
      write(45,'(A)')'DATASET STRUCTURED_GRID '
      write(45,'(A, I5, A1, I5, A1, I5)')
     & 'DIMENSIONS', dimZ+1, ' ', maxT+1, ' ',dimR+1
      write(45,'(A, I8, A)')'POINTS ',(dimZ+1)*(maxT+1)
     % *(dimR+1),' float'

      do iR = 0, dimR
        do iT = 0, maxT
          do iz = 0, dimZ

      posx = sin(dfloat(iT)/maxT*2.0*3.14159)*(iR)*delta
      posy = cos(dfloat(iT)/maxT*2.0*3.14159)*(iR)*delta
      posz = iz*delta

            write(45, *)
     & posx,'   ', ! sistema de coord x y
     & posy
     &, '   ', posz ! grafico en el viejo sistema x,y no en v,u
          enddo
        enddo
      enddo

      write(45,'(A, I8)')'CELL_DATA ', dimR*dimZ*maxT
      tempc = 'SCALARS ' // title // ' float 1'
      write(45,'(A)')tempc
      write(45,'(A)')'LOOKUP_TABLE default'

       do iR = 1, dimR
        do iT = 1, maxT
          do iZ = 1, dimZ

             single = array(iR, iZ) ! Lo necesito en single presicion

            write(45,*)single
          enddo
        enddo
      enddo
      close (45)

      endif

      if((savetodisk_type.eq.1).or.(savetodisk_type.eq.2)) then

      write(filename,'(A5, A1,I3.3, A1, I3.3, A4)')title,'.',
     &  counter, '.',  counter2, '.dat'

      open(unit=45, file=filename)

      do iR=1,dimR
         do iZ=1,dimZ
            write(45,*)(iR-0.5)*delta, (iZ-0.5)*delta, array(iR, iZ)
         enddo
      enddo

      close(45)

      endif

      return
      end


