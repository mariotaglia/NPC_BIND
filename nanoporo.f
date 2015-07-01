      program nanopore

      use mncells
      use mparameters
      use mparameters_chain
      use mparameters_monomer
      use mKaps

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      integer counter
      counter = 1

c--------------------------------------------------
!
!  Inits MPI
!

! MPI

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

       if(rank.eq.0)print*, 'Nuclear Pore Complex Program'
       if(rank.eq.0)print*, 'GIT Version: ', _VERSION
c--------------------------------------------------

       call globals       ! read global variables from file and allocate arrays
       if(rank.eq.0)print*, 'Globals OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call monomer_definitions ! set different types of monomer
      if(rank.eq.0)print*, 'Monomer definitions OK'
       if(rank.eq.0)print*, ' '
      call chains_definitions   ! set different types of chains and load protein sequences
      if(rank.eq.0)print*, 'Chains definitions OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call array_alloc ! array memory allocation
      if(rank.eq.0)print*, 'Memoray allocation OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call kai                 ! calculate or load kai
      if(rank.eq.0)print*, 'Kai OK'
       if(rank.eq.0)print*, ' '

      if(rank.eq.0) then
      call kap                 ! calculate or load kaps
      call kap_s               ! calculate or load kaps
      if(rank.eq.0)print*, 'Kaps OK'
       if(rank.eq.0)print*, ' '
      endif

      call geom                ! define pore shape
      if(rank.eq.0)print*, 'Geom OK'
       if(rank.eq.0)print*, ' '


      call lookup              ! define conectivity in matriz
      if(rank.eq.0)print*, 'Lookup OK'
       if(rank.eq.0)print*, ' '

      call create_protein
      if(rank.eq.0)print*, 'Create OK'
       if(rank.eq.0)print*, ' '


      CALL MPI_BARRIER(MPI_COMM_WORLD,err)
      call graftpoints         ! define graftpoint positions
      CALL MPI_BARRIER(MPI_COMM_WORLD,err)

      if(rank.eq.0)print*, 'Graftpoints OK'
       if(rank.eq.0)print*, ' '

      call creador             ! creates chains
      if(rank.eq.0)print*, ' '
      if(rank.eq.0)print*, 'Creador OK'

      call lookup_kai          ! sets up kai matrixes in lattice
      if(rank.eq.0)print*, 'Lookup_kai OK'
       if(rank.eq.0)print*, ' '

      if(rank.eq.0) then
      call lookup_kap          ! sets up kai matrixes in lattice
      call lookup_kap_s          ! sets up kai matrixes in lattice
      if(rank.eq.0)print*, 'Lookup_kap OK'
       if(rank.eq.0)print*, ' '
      endif

      call mainsolver(counter)          ! perform calculation and save results to disk

      call MPI_FINALIZE(ierr) ! finalize MPI

      end

      subroutine mainsolver(counter)

      use mlookup
      use mparameters ! parametros     
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mparameters_monomer
      use mprotein
      use mrands
      use mKaps

      implicit none

      include 'mpif.h'
      include 'MPI.h' ! MPI libraries

      integer i,j, k, ic, ii,ir,iz, im
      integer counter
      real*8 stok

      real*8 temp
      real*8 temp1

      character basura
      character*24 filename 
      character*5  title

      real*8 error
      real*8 sumpol, fmedio     

C-----  solving variables  -----------

      real*8 algo

      integer flag
       
      real*8 errel, fnorm
      integer n, itmax

      integer cc,ccc,cccc


C--------------------------------------------

! Kinsol

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations

! IMSL

      external fcnelect

C---------------------------------------------
      common /psize/ neq ! Kinsol
      external fcn


      real*8 x1((2+N_poorsol)*ncells)
      real*8 xg1((2+N_poorsol)*ncells)
      real*8 xflag((2+N_poorsol)*ncells)
      real*8 xflag2((2+N_poorsol)*ncells)
      real*8 f((2+N_poorsol)*ncells)

! auxiliary bulk

      real*8 prokap

! Variables

      shift = 1d0

      lb = 0.714 ! bjerrum lenght in nm

      zpos = 1.0
      zneg = -1.0
      
      vsol = 0.030     
      vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
      vpol= 0.060/vsol! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 

      constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  

      pKw = 14
      Kw = 1.0e-14

      error = 1e-6 ! para comparar con la norma...

      betae = 38.94             ! beta * e


      n=ncells
      
      errel=1d-6
      itmax=200

C-------------------------------------------------------
C MAIN LOOP
C-------------------------------------------------------

      flag = 0
      ccc = 1
      cc = 1

      do while(ccc.le.nxkapbulk) ! loop in st

 555  if(rank.eq.0)print*, 'ccc', ccc, 'of', nxkapbulk

       xkapbulk = xkapbulks(ccc)


c-------------------------------------------------------
c
c  CASE DEPENDENT VARIABLES
c

 257   cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
       xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
       pOHbulk= pKw -pHbulk
       cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
       xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  

! fdis bulk for the kap
       do im = 1, N_monomer
       Ka(im)=10**(-pKa(im))
       Kbind(im)=10**(-pKbind(im))
       Kbind0(im) = Kbind(im)*(1.0d24/Na)/vsol ! see notes
       enddo
   
       do im =1,N_monomer
           if (zpol(im).eq.1) then !BASE
              fdisbulk(im) =
     &        1.0 /(1.0 + xOHminbulk/((Kw/Ka(im))*vsol*Na/1.0d24))
            else if (zpol(im).eq.-1) then !ACID
               fdisbulk(im) =
     &         1.0 /(1.0 + xHplusbulk/(Ka(im)*vsol*Na/1.0d24))
            endif
       enddo

      vkap = sum(nmonkap)*vpol ! total volume of the kap in units of vsol 
      xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

      if(pHbulk.le.7) then  ! pH<= 7
            xposbulk=xsalt/zpos
            xnegbulk=
     &  - xsalt/zneg + (xHplusbulk-xOHminbulk)*vsalt ! NaCl+ HCl  
      else                  ! pH >7 
            xposbulk=xsalt/zpos +(xOHminbulk-xHplusbulk)*vsalt ! NaCl+ NaOH   
            xnegbulk=-xsalt/zneg 
      endif

      do im =1,N_monomer
           if (zpol(im).eq.1) then !BASE
              xnegbulk = xnegbulk 
     &      - xkapbulk/vkap*nmonkap(im)*fdisbulk(im)*vsalt*zpol(im)/zneg
            else if (zpol(im).eq.-1) then !ACID
              xnegbulk = xnegbulk 
     &      - xkapbulk/vkap*nmonkap(im)*fdisbulk(im)*vsalt*zpol(im)/zneg
            endif
      enddo

! xsol bulk

       xsolbulk=1.0 -xHplusbulk -xOHminbulk -
     &        xnegbulk -xposbulk - xkapbulk

! cargo pKa
         KaP=10**(-pKaP)
         if (prot_q < 0) then
         K0P = (KaP*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
         else ! base
         K0P = ((Kw/KaP)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
         endif

! Kaps K0
         do im = 1, N_monomer
         select case (zpol(im))
         case (-1) ! acid
         K0(im) = (Ka(im)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
         case (1) ! base
         K0(im) = ((Kw/Ka(im))*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
         end select
         enddo

! calculation xtotalbulk

         xtotalbulk = 0.0

         do im = 1, N_monomer
         ii = hydroph(im)
         xtotalbulk(ii) = xtotalbulk(ii)
     &   + xkapbulk/sum(nmonkap)*nmonkap(im)
         enddo

! calculation expmubulk

! 1. xpot

      do im =1, N_monomer
      xpotbulk(im) = dlog(xsolbulk)*vpol

      if(zpol(im).ne.0) then
       xpotbulk(im) = xpotbulk(im)
     & -dlog(fdisbulk(im))
      endif

      if(hydroph(im).ne.0) then
      do ii = 1, N_poorsol ! loop over different poor sv types
      xpotbulk(im) = xpotbulk(im) +
     &   (st_matrix(hydroph(im),ii) ! st_matrix(x, y) : interaction of hydrophobic segments of type x with those of type y 
     &   *st/(vsol*vpol)*           ! st in kT/monomer          
     &   sumXu*
     &   xtotalbulk(ii))
      enddo ! ii
      endif ! hydrophob
      enddo ! im

! 2. prokap
 
       prokap = 0.0

       do im = 1, N_monomer
           prokap = prokap 
     &   + xpotbulk(im)*nmonkap(im)
       enddo 

       expmukap = xkapbulk/dexp(prokap)

! chemical potentials

       expmupos = xposbulk /xsolbulk**vsalt
       expmuneg = xnegbulk /xsolbulk**vsalt
       expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
       expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

        if(rank.eq.0) then
         print*, 'Bulk composition / volume fracion'
         print*, 'Cations ', xposbulk
         print*, 'Anions  ', xnegbulk
         print*, 'H+      ', xHplusbulk
         print*, 'OH-     ', xOHminbulk
         print*, 'Solvent ', xsolbulk
         print*, 'Kaps    ', xkapbulk
         print*, 'Charge / vsol units'
 
         print*, 'Cations', xposbulk/vsalt*zpos         
         print*, 'Anions', xnegbulk/vsalt*zneg         
         print*, 'H+', xHplusbulk         
         print*, 'OH-', -xOHminbulk         

         temp = 0.0
         do im = 1, N_monomer
         temp = temp + 
     & xkapbulk/vkap*nmonkap(im)*fdisbulk(im)*zpol(im)
         print*, 'Kaps, type ',im, 
     & xkapbulk/vkap*nmonkap(im)*fdisbulk(im)*zpol(im)         
         enddo
 
         print*, 'Kaps', temp

         temp = temp +        
     & xposbulk/vsalt*zpos+xnegbulk/vsalt*zneg+xHplusbulk-xOHminbulk

         print*, 'all', temp

         print*, 'Monomer type', ' K0 ','fdisbulk'
         do im = 1, N_monomer
         if(zpol(im).ne.0)print*, im, K0(im), fdisbulk(im) 
         enddo

         print*, 'Poor sv type ', ' xtotalbulk'
         do ii = 1, N_poorsol
         print*, ii, xtotalbulk(ii)
         enddo

         print*, 'vkap (in nm^3)', vkap*vsol
         print*, 'Kap volume from diameter',
     &   (4.0/3.0)*pi*((float(Kapd)/2.0)**3)*delta**3

         endif

!         call MPI_FINALIZE(ierr) ! finalize MPI
!         stop

c-------------------------------------------------------------------------------
c
c INITIAL GUESS 
c
c-------------------------------------------------------------------------------

      if(infile.eq.2) then
      if(ccc.ne.1) then
      do i = 1, (2+N_poorsol)*ncells ! dejar aca, asegura guess correcto para infile = 2 (cc > 1)
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
      enddo
      endif
      endif
c
      if(infile.eq.0) then
      do i=1,n
         xg1(i)=xsolbulk
         x1(i)=xsolbulk
      enddo

      do i=n+1, n*2
         xg1(i)=0.0d0
         x1(i)=0.0d0
      enddo

      do i=(1+N_poorsol)*ncells+1, n*(2+N_poorsol)
          xg1(i)=0.0
          x1(i)=0.0
      enddo

      endif
 
      if(infile.eq.1) then
      open(unit=45, file='in.txt')
      do i=1, (2+N_poorsol)*ncells
      read(45,*)xg1(i)
      x1(i) = xg1(i)
      enddo
      close(45)
      endif 

      if(infile.eq.5) then
      write(filename,'(A4, I3.3, A4)')'out.', ccc, '.dat'
      print*, 'reading ', filename
      open(unit=45, file=filename)
      do i = 1, (2+N_poorsol)*ncells
      read(45, *)x1(i)
      enddo
      xg1 = x1
      close(45)
      endif

C--------------------------------------------------------------
C               +++ SOLVER +++ 
C--------------------------------------------------------------

! JEFE

          if(rank.eq.0) then ! solo el jefe llama al solver

          iter = 0
          print*, 'Enter Solver ', ncells*(2+N_poorsol), ' equations.'
            

          if(infile.ne.5) then       
          call call_kinsol(x1, xg1, ier)
          else
          call fkfun(x1, f, ier) ! data analisis
          endif
          
          flagsolver = 0
          
          CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, 
     & MPI_COMM_WORLD,err)
          endif

! Subordinados

      if(rank.ne.0) then

      do

        flagsolver = 0
        source = 0

       CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)

        if(flagsolver.eq.1) then
       CALL MPI_BCAST(x1, (2+N_poorsol)*ncells , MPI_DOUBLE_PRECISION,0
     &  , MPI_COMM_WORLD,err)

        call fkfun(x1, f, ier) ! todavia no hay solucion => fkfun 
        endif ! flagsolver

        if(flagsolver.eq.0)goto 1010 ! Detiene el programa para este nodo

        enddo

        endif

 1010   continue
          
! Retrives ier and nomr
! This way slaves know if solver process in master has converged
!

          ! master

          if (rank.eq.0) then

            norma_tosend = norma
            ier_tosend = ier 

            CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,  
     &      MPI_COMM_WORLD,err)

            CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER, 0,
     &      MPI_COMM_WORLD,err)

         endif

         ! slaves

         if (rank.ne.0) then

            CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,  
     &      MPI_COMM_WORLD,err)

            CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER, 0,
     &      MPI_COMM_WORLD,err)

         norma = norma_tosend
         ier = ier_tosend

         endif
          
! Retrives xsol y psi (NOT COMMON!)

            do iC = 1, n
            xsol(iC)=x1(iC)
            psi2(iC)=x1(iC+n)
            do i = 1, N_poorsol
            xtotal2(i, iC)=x1(iC+(i+1)*n)
            enddo
            enddo

            xtotal2(:,ncells+1)=xtotalbulk(:)
            xtotal2(:,ncells+2)=xtotalbulk(:)

! End of solver

! Checks for convergence

         if(infile.ne.5) then
         if((ier.lt.0).or.
c     &   (not((norma.gt.0).or.(norma.lt.0))).or. ! remove this line when compiling with gfortran
     &   (norma.gt.error)) then ! not converged

         if(rank.eq.0) then
         print*, 'Error in the solver: ', ier
         print*, 'norm ', norma
         print*, 'Error xkapbulk', xkapbulk
         endif
        
         if(ccc.ne.1) then
         algo = (xkapbulk + stok)/2
         else
         algo = xkapbulk/2
         endif

c         nsigma = nsigma + 1

         if(rank.eq.0) then
         write(1010,*)'xkapbulk failed ', xkapbulk, ' go to ', algo
         endif

         xkapbulk = algo
         flag = 1

         goto 257
         endif    
         endif

! Solution converged, store xflag

         do i = 1, (2+N_poorsol)*n
         xflag(i) = x1(i) ! xflag will be used as input for the next iteration
         end do

         if(infile.ne.5)infile = 2 ! no vuelve a leer infile
         stok = xkapbulk

         if(flag.eq.1) then  ! recorvers from an error
         if(rank.eq.0) then
         print*, 'Recovers from an error'
         print*, 'OK', stok
         endif
         flag = 0
         goto 555
         endif



      if(rank.eq.0) then ! only master can save to disk

! saves infile

      if(infile.ne.5) then
      write(filename,'(A4, I3.3, A4)')'out.', ccc, '.dat'
      open(unit=45, file=filename)
      do i = 1, (2+N_poorsol)*ncells
      write(45, *)x1(i)
      enddo
      close(45)
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cc = saveindex 
      call saveresults(ccc, cc) ! save results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call calc_free_energy(dfloat(ccc), dfloat(cc)) ! calculate free energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ccc = ccc + 1

      enddo ! ccc

      end


