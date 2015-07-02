c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)

      use mparameters
      use mkinsol
      use mparameters_chain
      use mparameters_monomer

      implicit none

      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)

      common /psize/ neq

      do  i = 1, neq
         vv(i) = vv(i) * pp(i)
      enddo
      ier = 0

      return
      end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

      use mparameters
      use mkinsol
      use mparameters_chain
      use mparameters_monomer

      implicit none

      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      integer n

      common /psize/ neq

      n = neq/(2+N_poorsol)

      do i = 1, n
         pp(i) = 1.0 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
      enddo

      do i = n+1,2*n
         pp(i) = 1.0
      enddo

      do i = 2*n+1, neq
         pp(i) = -1.0 / (1.0+exp(-udata(i)))
      enddo



      ier = 0

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subrutina que llama a kinsol
!

      subroutine call_kinsol(x1, xg1, ier)

      use mncells
      use mparameters
      use mparameters_chain
      use mparameters_monomer

      implicit none

      real*8 x1((2+N_poorsol)*ncells), xg1((2+N_poorsol)*ncells)

      integer*8 iout(15) ! Kinsol additional output information
      real*8 rout(2) ! Kinsol additional out information

      integer*8 msbpre
      real*8 fnormtol, scsteptol
      real*8 scale((2+N_poorsol)*ncells)
      real*8 constr((2+N_poorsol)*ncells)

      integer*4  globalstrat, maxl, maxlrst, max_niter

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations

      integer i, ierr

      common /psize/ neq ! Kinsol


! INICIA KINSOL


      neq = (2+N_poorsol)*ncells

     

      msbpre  = 1000 ! maximum number of iterations without prec. setup (?)
      fnormtol = 1.0d-8 ! Function-norm stopping tolerance
      scsteptol = 1.0d-8 ! Function-norm stopping tolerance

      maxl = 5000 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
      maxlrst = 200 ! maximum number of restarts

      max_niter = 2000

      globalstrat = 1

      call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
      if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
      print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
      call MPI_FINALIZE(ierr) ! finaliza MPI
         stop

      endif

      call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
      if (ier .ne. 0) then
         print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
         call MPI_FINALIZE(ierr) ! finaliza MPI
         stop
      endif


      call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      call fkinsetrin('SSTEP_TOL', scsteptol, ier)
      call fkinsetiin('MAX_NITER', max_niter, ier)
                                                      
         do i = 1, ncells  !constraint vector
         constr(i) = 2.0
         enddo
         do i = ncells+1, 2*ncells  !constraint vector
         constr(i) = 0.0
         enddo
         do i = 2*ncells+1, neq
         constr(i) = 0.0
         enddo 

         call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector

c       CALL FKINSPTFQMR (MAXL, IER)
        call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
c       call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG


      if (ier .ne. 0) then
      print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
         call fkinfree ! libera memoria
         call MPI_FINALIZE(ierr) ! finaliza MPI
         stop
      endif
      call fkinspilssetprec(1, ier) ! preconditiones

      do i = 1, neq ! scaling vector
      scale(i) = 1.0
      enddo

    

      do i = 1, neq ! Initial guess
      x1(i) = xg1(i)
      enddo


      call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol
      if (ier .lt. 0) then
      print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'Linear Solver returned IER = ', iout(9)
c         call fkinfree

      endif

      return
      end



