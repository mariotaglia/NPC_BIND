      subroutine creador ! crea cadenaes

      use mlookup
      use mchains
      use mparameters
      use mncells
      use mgraft
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mprotein

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
 
      integer k,vx(4),vy(4)

      integer total,ix(3)
      integer *1 displ_temp
      integer i,il,u1,u2,iii,ii,ll, jj
      integer j,ncha

      real*8 indax, inday, indaz
      
      real*8 altx,alty, altz
      real*8 rij,theta,theta1, rn1, rn2

      real*8 rpos, zpos
           
      integer total1,iglobal

      integer spp ! sitios por procesador

      integer cuantas_p(2)

      real*8 q_tosend
      integer*1 displ_one(0:7)

      integer*1 displacement(maxlong+1,2)
      integer*1 binary(int(maxlong/2))
     
      real*8 xend(3,maxlong)
      integer chains(ncha_max, maxlong)
      real*8 endtoendtemp(10000)
      real*8 x(maxlong)
      real*8 y(maxlong)
      real*8 z(maxlong)

      real*8 in1(maxlong, 3)
      common /endtoend/ endtoendtemp

      seed = readseed
 
      rpos = posgraft(1)
      zpos = posgraft(2)

      ii = rank+1
    
      newcuantas(ii) = 0

      il=0
      iglobal=1
 
      if(calq.eq.1)ncha=-1

      do while (il.lt.cuantas)
         
         if(cadenastype.eq.1)call cadenas72mr(chains,ncha, rpos, zpos, 
     & long(ii))
         if(cadenastype.eq.2)
     & call cadenas_mk(chains,ncha, rpos, zpos, long(ii))
      
         do i=1,ncha

         il = il + 1
 
         if(il.gt.cuantas) goto 100

! Check collision with protein...

         do j = 1, long(ii)
         if(proteinC(chains(i,j)).ne.0.0)goto 200 ! collides with protein, don't use this chain
         enddo

         newcuantas(ii) = newcuantas(ii)+1

         endtoend(newcuantas(ii)) = endtoendtemp(i)

!          fs(newcuantas(ii))=chains(i,1)

!          if(fs(newcuantas(ii)).eq.0) then 
!          print*, 'Error', il, i, rank, ncha
!          stop
!          endif
 
           do j = 1, long(ii)
           inc(newcuantas(ii),j) = chains(i,j)
           enddo
         

!          displacement(1,1) = 0
!          displacement(1,2) = 0

!          do j=2,long(ii)        
!
!          displacement(j,1)= 
!     &  (indexa(chains(i,j),1)-indexa(chains(i,j-1),1)) ! displacement in R
          
!          displacement(j, 2) =
!     &  (indexa(chains(i,j),2)-indexa(chains(i,j-1),2)) ! displacement in Z


!          if(displacement(j,2).eq.(1-dimZ))displacement(j,2)=1    ! PBC
!          if(displacement(j,2).eq.(-1+dimZ))displacement(j,2)=-1

!          enddo ! 

!          call encode(displacement,binary, long(ii))
!          displ(newcuantas(ii),:) = binary(:)

 200      enddo !ncha
          enddo ! il 

      print*, 'Processor ', rank+1, 'has',newcuantas(ii),'conformations'
      if(newcuantas(ii).eq.0)stop
          
 100  return
      end

C****************************************************************

      subroutine cadenas72mr(chains,ncha, rpos, zpos, long1)

      use mchains
      use mparameters
      use mncells
      use msegme
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mporesystem

      implicit none

      integer ncha
      integer chains(ncha_max,maxlong)
      integer long1
      
      real*8 y(maxlong),z(maxlong)
      
      real*8 rvect, zvect
      
      integer i,state,ii,j,ive,jve
      real*8 rn,state1,sitheta,cotheta,dista

      real*8 siphip,cophip
      character*1 test
      real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
      
      real*8 x(3),xend(3,maxlong),xendr(3,maxlong), xendr2(3,maxlong)
      real*8 theta
      real*8 zpos, rpos
      
      real*8 tempp(3)

      logical   outside
      external  outside

      integer iR, iZ

      real*8 xo(3)
      
      sitheta=sin(68.0*pi/180.0)
      cotheta=cos(68.0*pi/180.0)
      siphip=sin(120.0*pi/180.0)
      cophip=cos(120.0*pi/180.0)

      theta = pi
            
 223  x(1)=lseg
      x(2)=0.0
      x(3)=0.0
      
      xend(1,1)=lseg
      xend(2,1)=0.0
      xend(3,1)=0.0
      
      tt(1,1)=cotheta
      tt(1,2)=sitheta
      tt(1,3)=0.0
      tt(2,1)=sitheta
      tt(2,2)=-cotheta
      tt(2,3)=0.0
      tt(3,1)=0.0
      tt(3,2)=0.0
      tt(3,3)=-1.0
      
      tp(1,1)=cotheta
      tp(1,2)=sitheta
      tp(1,3)=0.0
      tp(2,1)=sitheta*cophip
      tp(2,2)=-cotheta*cophip
      tp(2,3)=siphip
      tp(3,1)=sitheta*siphip
      tp(3,2)=-cotheta*siphip
      tp(3,3)=-cophip
      
      tm(1,1)=cotheta
      tm(1,2)=sitheta
      tm(1,3)=0.0
      tm(2,1)=sitheta*cophip
      tm(2,2)=-cotheta*cophip
      tm(2,3)=-siphip
      tm(3,1)=-sitheta*siphip
      tm(3,2)=cotheta*siphip
      tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
      state1=0.0
      
      m(1,1)=cotheta
      m(1,2)=sitheta
      m(1,3)=0.0
      
      m(2,1)=cos(state1)*sitheta
      m(2,2)=-cos(state1)*cotheta
      m(2,3)=sin(state1)
      m(3,1)=sin(state1)*sitheta
      m(3,2)=-sin(state1)*cotheta
      m(3,3)=-cos(state1)
      
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
      
      xend(1,2)=lseg+x(1)
      xend(2,2)=x(2)
      xend(3,2)=x(3)

      do 10 i=3,long1
         
 123     rn=rands(seed)
         state=int(rn*3)
c     print*,'state',state

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
            
            call mrrrr(m,tt,mm)
            do 30 ii=1,3
               do 40 j=1,3
                  m(ii,j)=mm(ii,j)
 40            continue
 30         continue
            
            
         elseif (state.eq.1) then
            
            call mrrrr(m,tp,mm)
            do 31 ii=1,3
               do 41 j=1,3
                  m(ii,j)=mm(ii,j)
 41            continue
 31         continue

         elseif (state.eq.2) then
            
            call mrrrr(m,tm,mm)
            do 32 ii=1,3
               do 42 j=1,3
                  m(ii,j)=mm(ii,j)
 42            continue
 32         continue
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
         
c     if (xend(1,i).lt.0.0) goto 222
         
 10   continue
      
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chequea cadenas
! Selfavoiding entre segmentos
 
      dista=0.0
      do 300 ive=4,long1
         do 310 jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
            
               goto 222
            endif
 310     continue
 300  continue

      ncha=0
      do ii=1,300

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random rotation
!

      call rota36(xend,xendr,long1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Displace chains to (rpos, zpos)
! 
        
      tempp(1) = xendr(1, 1)
      tempp(2) = xendr(2, 1)
      tempp(3) = xendr(3, 1)
 
      do i=1,long1 
 
         xendr(1, i) =  xendr(1, i) + rpos - tempp(1) - 1e-4
         xendr(2, i) = xendr(2, i) - tempp(2)
         xendr(3, i) =  xendr(3, i) + zpos - tempp(3)
       
         
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
! Put segment in lattice / check coordintes      
!     
         do i = 1, long1
         do j = 1, 3
         in1(i, j) = xendr(j, i)
         xo(j) = xendr(j, i)
         enddo
         if(outside(xo)) goto 400
         enddo

         ncha=ncha+1

         do j = 1, long1

         do while(in1(j,3).ge.dfloat(dimZ)*delta)
           in1(j,3) = in1(j,3) - dfloat(dimZ)*delta
         enddo
          do while(in1(j,3).le.0.0)
           in1(j,3) = in1(j,3) + dfloat(dimZ)*delta
         enddo

         iR = int(sqrt(in1(j,1)**2+in1(j,2)**2)/delta)+1
         iZ = int(in1(j,3)/delta)+1

         if (iR.gt.dimR) then
         print*, 'Error in creador'
         print*, 'Increase system size'
         print*, iR, iZ
         stop
         endif

         chains(ncha,j) = matriz(iR,iZ)
  
         enddo
     
         if (ncha.ge.25) goto 402
         
 400  enddo
 402  if (ncha.eq.0) goto 223

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rota36(xend,xendr,n)
      
      use  mparameters
      use mparameters_chain
      use mparameters_monomer
      use mrands

      implicit none
 
      real*8 xend(3,maxlong),xendr(3,maxlong)
      character*1 test
      integer n
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      
      real*8 radio, vect
      real*8 alfa, gama, cga, a, b, c
      integer i
 

      
      fac=rands(seed)
      fac1=rands(seed)
      fac2=rands(seed)
      alfa=fac*2*pi
      cbe=2.0d0*fac1-1.0d0
      gama=fac2*2*pi

      sbe=(1-cbe**2)**0.5
      cal=cos(alfa)
      sal=sin(alfa)
      cga=cos(gama)
      sga=sin(gama)

      do 1 i=1,n

         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)
     &        -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+
     &        b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
 
 1    continue

      return
      end
      
C****************************************************************

      subroutine mrrrr(a,b,c)

      implicit none
      
      real*8 a(3,3),b(3,3),c(3,3)

      integer i, j, k

      do 1 i=1,3
         do 1 j=1,3
            c(i,j)=0
 1    continue

      do 2 i=1,3
         do 2 j=1,3
            do 2 k=1,3
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2    continue

      return
      end 

      subroutine cadenas_mk(chains,ncha,rpos,zpos,N)
! =====================================================================================================
! version 21 may 2011 mk@mat.ethz.ch
! =====================================================================================================

! =====================================================================================================
! cadenas requires ncha = -1 before it is called for the first time! 
! =====================================================================================================

! creates ncha polymers of length N confinded via matriz >= 0
! returns bead positions as cell numbers (chains)
! bond length: lseg, bond diameter: dseg (the latter used in overlap check)
! chains(1..ncha,1..long) = matriz(iR,iZ)

! matriz(iR,iZ) = cell number (outside geometry: < 0)   ! pore.system.h
! indexa(cell number,1) = iR                            ! lookup.h
! indexa(cell number,2) = iZ

        use mparameters
        use mporesystem
        use mparameters_chain
        use mparameters_monomer

        implicit none
        include 'mpif.h' ! MPI libraries
        include 'MPI.h' ! MPI libraries

        real*8 rpos,zpos,x(3)
        real*8 endtoendtemp(10000)
        integer ncha,chains(ncha_max,maxlong),N
        integer ncha_current
        common /comncha/ ncha_current
        integer N_max
        common /compass/ N_max
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /endtoend/ endtoendtemp

        real*8 qprob
        common /qprob/ qprob

! =====================================================================================================
        branch_p      = 5               ! EVENTUALLY SYSTEM AND MODEL SPECIFIC
        spacer        = 2               ! BUT q IS AUTOMATICALLY ADJUSTED TO CHOSEN s and p VIA 'TUNING'
! =====================================================================================================

! =====================================================================================================
        x(1) = rpos
        x(2) = 0.D0
        x(3) = zpos
! =====================================================================================================
! =============================== DO NOT EDIT BELOW THIS LINE =========================================
! =====================================================================================================

        ncha_current    = 0
        N_max           = N
        if ((calq.eq.1).and.(ncha.eq.-1)) then 
        call TUNING(N,rpos,zpos)
        call MPI_FINALIZE(ierr) ! finalize MPI
        stop
        endif

        call walk(1,x,chains)
        ncha = ncha_current
        return
        end

! =========================================================================== 

        recursive subroutine walk(N,x,chains)
        use mparameters
        use mgraft
        use mporesystem
        use mrands
        use mparameters_chain
        use mparameters_monomer
        use posmk

        implicit none

        real*8 x(3),y(3),u(3)
        integer N
        integer ncha_current
        common /comncha/ ncha_current
        integer*2 beadno
        integer   j,k,ix1,ix2,ix3,iR,iZ,ibranch
        real*8    dummy,stretched,power
        logical   outside,hit_bead
        external  outside,hit_bead
        integer N_max
        common /compass/ N_max
        integer chains(ncha_max,maxlong)
        real*8 endtoendtemp(10000)
c        integer*4, allocatable ::  deathcount(:)
        integer*4 deathtotal
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /crit/ deathtotal
        real*8 qprob
        common /endtoend/ endtoendtemp
        common /qprob/ qprob
        save power
c        save deathcount

c        if(.not.allocated(deathcount))ALLOCATE(deathcount(maxlong)) 

        if (N.eq.1) then
         ! print 3,N,'root',ncha_current
         stretched   = (N_max-1)*lseg                           ! EVENTUALLY MODEL SPECIFIC
         power = dlog(dble(mcube))/dlog(stretched/dseg)
         power = min(1.D0,power)
         ncha_current = 0
         firstcell    = 0
         nextbead     = 0
         beadno       = 0
         nextbead(N)  = 0
         firstcell(0,0,0) = N
c         deathcount   = 0
         deathtotal   = 0
         ! print 1,'long',long
         ! print 1,'N',N
         ! print 2,'qprob0',qprob0
         ! print 1,'spacer',spacer
         ! print 1,'branch_p',branch_p
         ! print 1,'nearbonds',nearbonds
         ! print 1,'ncha_current',ncha_current
         do k=1,3; current(1,k) = x(k); enddo
         if (outside(x)) then
          print *,x
          print *,int(sqrt(x(1)**2+x(2)**2)/delta),int(x(3)/delta)
          stop 'anchor is outside!?'
         end if
111      continue
         call randomunit(u)
         do k=1,3; y(k) = x(k) + lseg*u(k); enddo
         if (outside(y)) goto 111
         call walk(2,y,chains)          ! calls N=2, next gen is 2 (after some strand ..)
         return
        end if

        if (calq.eq.0) then

         if (rands(seed).gt.qprob0) then
          ! print 3,N,'q prob return',ncha_current
          goto 222
          return
         endif
         endif

        if (calq.eq.1) then

         if (rands(seed).gt.qprob) then
          ! print 3,N,'q prob return',ncha_current
          goto 222
          return
         endif

        endif


        if (ncha_current.ge.ncha_max)   return 
        if (outside(x)) then
         ! print 3,N,'outside return',ncha_current
         goto 222
         return
        endif

        dummy = x(1)-current(1,1)
        ix1 = nint(dexp(power*
     >        dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
        dummy = x(2)-current(1,2)
        ix2 = nint(dexp(power*
     >        dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
        dummy = x(3)-current(1,3)
        ix3 = nint(dexp(power*
     >        dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))

c        print*, 'firstcell is', size(firstcell, 1), size(firstcell, 2),
c     & size(firstcell, 3)

        if (hit_bead(ix1,ix2,ix3,x,N)) then
         ! print 3,N,'hitbead return',ncha_current
         goto 222
        end if

        beadno                   = firstcell(ix1,ix2,ix3)
        nextbead(N)              = beadno
        firstcell(ix1,ix2,ix3)   = N

        do k=1,3; current(N,k) = x(k); enddo

        if (N.eq.N_max) then                   ! complete chain generated successfully

         ncha_current = ncha_current + 1
         ! print 3,N,'chain end ---',ncha_current

         do j=1,N
          iR = int(sqrt(current(j,1)**2+current(j,2)**2)/delta)+1
          iZ = int(current(j,3)/delta)+1
         if ((iR.gt.dimR).or.(iZ.lt.1).or.(iZ.gt.dimZ)) then
         print*, 'Error in creador'
         print*, 'Increase system size'
         endif
          chains(ncha_current,j) = matriz(iR,iZ)
         enddo

         endtoendtemp(ncha_current) = 
     & ((current(N,1)-current(1,1))**2 +
     & (current(N,2)-current(1,2))**2 +
     & (current(N,3)-current(1,3))**2)**(0.5)

         ! call checking_actual_config(N)       ! DEACTIVE IN PRODUCTION

         if (ncha_current.ge.ncha_max) return

        else                                    ! N < N_max

         if (mod(N-1,spacer).eq.1) then         ! N=2,2+spacer,...
          ! print 3,N,'branching',ncha_current
          do ibranch=1,branch_p
           call randomunit(u)
           do k=1,3; y(k) = x(k) + lseg*u(k); enddo
           call walk(N+1,y,chains)              
           if (ncha_current.ge.ncha_max) return
          enddo
         else
          ! print 3,N,'linear growth',ncha_current
          call randomunit(u)
          do k=1,3; y(k) = x(k) + lseg*u(k); enddo
          call walk(N+1,y,chains)     
          if (ncha_current.ge.ncha_max) return
         end if

        end if

        firstcell(ix1,ix2,ix3) = beadno         ! set free
        return

1       format("walk ",A30,I10)
2       format("walk ",A30,F10.3)
3       format(I5,1x,A20,I10)

222     continue
c        deathcount(N) = deathcount(N) + 1
        deathtotal    = deathtotal + 1

        return
        end


        function outside(x)
        use mparameters
        use mporesystem
        implicit none
        logical outside
        real*8 x(3),r, zz
        integer IR,IZ   
         outside = .true.

         zz=x(3)

         do while(zz.gt.dfloat(dimZ)*delta)
           zz = zz - dfloat(dimZ)*delta
         enddo
          do while(zz.le.0.0)
           zz = zz + dfloat(dimZ)*delta
         enddo

         r  = sqrt(x(1)**2+x(2)**2)
         IR = int(r/delta)+1
         IZ = int(zz/delta)+1

         if (matriz(IR,IZ).gt.0) outside = .false.      
        return
        end

        subroutine randomunit(seg)                      
        use mrands
        implicit none
        real*8 seg(3),znorm,z(2),mysqrt
         znorm=2.D0
         do while (znorm.ge.1.D0)
          z(1) = 1.D0-2.D0*rands(seed)
          z(2) = 1.D0-2.D0*rands(seed)
          znorm = z(1)*z(1)+z(2)*z(2)
         enddo
         mysqrt = dsqrt(1.D0-znorm)
         seg(1) = 2.D0*z(1)*mysqrt
         seg(2) = 2.D0*z(2)*mysqrt
         seg(3) = 1.D0-2.D0*znorm
        return
        end
        
        function hit_bead(ix1,ix2,ix3,x,N)
        use mparameters
        use mparameters_chain
        use mparameters_monomer
        use posmk 

        implicit none

        logical hit_bead
        integer ix1,ix2,ix3,k,N,k1,k2,k3
        real*8 x(3)
        integer*2 beadno
        real*8 isd2
         hit_bead = .false.
         do k1=-1,1
         do k2=-1,1
         do k3=-1,1
c         print*, ix1, ix2, ix3, size(firstcell, 1), 
c     & size(firstcell, 2), size(firstcell, 3)
          beadno = firstcell(ix1+k1,ix2+k2,ix3+k3)
          if (beadno.ge.N) stop 'BAD'
          do while (beadno.gt.0)
           isd2 = (current(beadno,1)-x(1))**2
     >           +(current(beadno,2)-x(2))**2
     >           +(current(beadno,3)-x(3))**2
           if (isd2.lt.d2) then                                  ! far away sphere have diameter d
            if (abs(beadno-N).le.1) then                         ! treat adjacent as nonoverlapping 23 MAY 2011
             hit_bead=.false.
            elseif (abs(beadno-N).le.nearbonds) then
             if (isd2.lt.b2) then                                ! near spheres have diameter b
              hit_bead=.true.;  return
             end if
            else
             hit_bead=.true.;  return
            end if
           end if
           beadno = nextbead(beadno)
          enddo
         enddo
         enddo
         enddo
        return
        end


! NOT IN USE 
        subroutine checking_actual_config(N)
        use mparameters
        use mgraft
        use mporesystem
        use mrands
        use mparameters_chain
        use mparameters_monomer
        use posmk

        implicit none
        integer ncha

        integer N,k
        integer ncha_current
        common /comncha/ ncha_current
        integer i,j1,j2
        real*8  dist2
        print *,'checking .. ',ncha_current
        do j1=1,N-1
        do j2=j1+1,N            
         dist2 = 0.D0
         do k=1,3
          dist2 = dist2 + (current(j1,k)-current(j2,k))**2
         enddo
         if (dist2.lt.d2) then
          print *,j1,j2
          print *,dist2,' < ',d2
          stop 'BUG'
         end if 
        enddo
        enddo
        return
        end


! ================================================================== MK SPECIAL REQUIRES s and p, finds q
        subroutine TUNING(N,rpos,zpos)

        use mlookup
        use mparameters
        use mporesystem
        use mrands
        use mparameters_chain
        use mparameters_monomer

        implicit none

        integer N,ncha
        real*8 tic,q_best
        integer cuantas_done,check_deaths,check_rounds
        real*8 cputotal,check_DAR,check_ok
        integer check_CPT
        real*8 rpos,zpos
        integer j
        integer chains(ncha_max,maxlong)
        integer*4 deathtotal
        common /crit/ deathtotal
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /qprob/ qprob

        real*8 qprob,deltaq
        logical converged
        print '(A,I6,A)','+++ TUNING wanted ',wantedCPUsecs,' secs'
        ! print *,'+++ mk-tuning cpu ',wantedCPUsecs/dble(1000)
c        print 101,'N','s','p','q','DAR','CPT','cpu','cpu-total','ok%'
        qprob  = 1.0            ! DO NOT EDIT
        deltaq = 0.1            ! DO NOT EDIT
        converged = .false.
        ncha   = 1
        do while (.not.converged)
         tic = secnds(0.)
         cuantas_done = 0
         check_deaths = 0
         check_rounds = 0
         do while ((secnds(0.).lt.tic+wantedCPUsecs/1000).      ! DO NOT EDIT
     >             or.(check_rounds.lt.40))                     ! DO NOT EDIT
          call cadenas_mk(chains,ncha,rpos,zpos,N)
          check_rounds = check_rounds + 1
          cuantas_done = cuantas_done + ncha
          check_deaths = check_deaths + deathtotal
         enddo
         tic = secnds(0.)-tic
         cputotal = tic*cuantas/dble(1e-6+cuantas_done)
         tic = 10**6*tic/dble(max(1,cuantas_done))
         check_DAR = check_deaths/dble(max(1,cuantas_done))
         check_CPT = cuantas_done/dble(check_rounds)
         check_ok  = 100*abs(cputotal-wantedCPUsecs)/dble(wantedCPUsecs)
c         print 100,N,spacer,branch_p,qprob,check_DAR,
c     >         check_CPT,tic,cputotal,check_ok
         if (int(cputotal).gt.wantedCPUsecs) goto 1
         if (cuantas_done.eq.0) goto 1
         qprob = qprob - deltaq
         goto 2
1        continue
         qprob = min(1.D0,qprob + deltaq)
         deltaq = deltaq/1.5
         if (deltaq.lt.1e-4) goto 3
2        continue
        enddo
3       continue
        print *,'now operating at N, q ',N,qprob
        

c100     format("+++ ",3(I6,1x),F8.6,1x,F8.3,1x,I8,1x,F8.1,2(1x,F8.1))
c101     format("+++ ",3(A6,1x),A8,1x,A8,1x,A8,1x,A8,1x,A8,1x,A8)

        return
        end


