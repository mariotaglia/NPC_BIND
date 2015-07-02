      !
      ! This rountine performs integration of an sphere in order to determine the contribution of its volume to different layers
      ! The result is stored in the array Ka 




      subroutine kap

      use mvariables
      use mparameters
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mKaps
  

      implicit none

      include 'mpif.h'
      include 'MPI.h'
                                                     ! interaction of a
                                                     ! segment at R = a
                                                     ! and z = 0, with a
                                                     ! segment at R = b
                                                     ! and Z = c
      real*8 R,theta,z
      real*8 thetamin, thetamax, rmax, rmin, zmax, zmin, r0
      real*8 cutoff

      real*8 rn
      integer i, ii
      real*8 radio ! total radius
      real*8 x1,x2,y1, y2, z1, z2, vect
      integer Rj, Zj, j
      character*40 filename

      real*8 rounding

      integer im
      real*8 vsphere

      integer iR,iZ,itheta

      real*8 sumtemp

      vsphere = 4.0/3.0*pi*(Kapd*delta/2.0)**3
      if (readkap.ne.1) then ! do not read from file

         radio = delta*dimR
         if(rank.eq.0)
     & print*,'Sphere calculation, readkap =', readkap

         write(filename,'(A5, I4.4,A1,I2.2, A4)')
     & 'Kaps-', dimR,'-',Kapd, '.kai'

         open(unit=112, file=filename)

         KapsV = 0.0               ! vector Kaps
         seed = readseed

         cutoff = Kapd*delta/2.0

         rounding = 1e-6 ! used to prevent errors due to rounding

         do ii = 1, dimR ! center of the particle    

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

!            if(rank.eq.0) then
!            print*, ii
!            print*, rmin, rmax
!            print*, thetamin, thetamax
!            print*, r0
!            endif          
!            thetamin = 0.0
!            thetamax = 2.0*pi           


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

               if((Zj.ge.((Kapd-1)/2+1)).or.(Zj.le.(-(Kapd-1)/2-1)))then
                print*, z, z/delta,i,ii
                stop
               endif

               if(vect.gt.(cutoff))cycle ! outside sphere

               if(Rj.gt.dimR)Rj=dimR+1 ! goes to bulk

               KapsV(ii, Rj, Zj) = KapsV(ii, Rj, Zj) + R ! incluye el jacobiano R(segmento)

            enddo ! iR
            enddo ! iZ
            enddo ! itheta

            sumtemp = 0.0
            do Rj = 1, dimR+1
            do Zj = -((Kapd-1)/2),((Kapd-1)/2)  ! Maximum z = (Kapd-1)/2
            sumtemp = sumtemp + KapsV(ii,Rj,Zj)
            enddo 
            enddo

            do Rj = 1, dimR+1
            do Zj = -((Kapd-1)/2),((Kapd-1)/2)  ! Maximum z = (Kapd-1)/2
            KapsV(ii, Rj, Zj) = KapsV(ii, Rj, Zj)/sumtemp*vsphere
!            KapsV(ii, Rj, Zj)=KapsV(ii, Rj,Zj)/(MCsteps**3)*(zmax-zmin)*
!     &                 (thetamax-thetamin)*(rmax-rmin)
            write(112,*)ii,Rj,Zj,KapsV(ii,Rj,Zj)
            enddo
            enddo

         end do                 ! ii -> Particle position
         close(112)
      else                      ! 
         write(filename,'(A11, I4.4,A1,I2.2, A4)')
     & '../../Kaps-', dimR,'-',Kapd, '.kai'
         open(unit=112, file=filename)
         if(rank.eq.0)print*,'Read Spheres'
         do j = 1, dimR*(dimR+1)*Kapd
            read(112,*)ii,Rj,Zj,KapsV(ii,Rj,Zj)

! OJO
!         if((ii.eq.Rj).and.(Zj.eq.0)) then
!         KapsV(ii,Rj,Zj) = 4.0/3.0*pi*(Kapd*delta/2.0)**3 
!         else
!         KapsV(ii,Rj,Zj) = 0.0
!         endif

         enddo
      endif
 
!! Assign segments according to the volume distribution of the Kap
!! Right now the volume assignement is done by distributing all types of different monomers 
!! proportional to the volume of the sphere 


      do im = 1, N_monomer
      if(center(im).ne.0) then
      Kaps(:,:,:,im) = KapsV(:,:,:)*nmonkap(im)/vsphere
      else
      Kaps(:,:,:,im) = 0.0
      do ii = 1, dimR
      Kaps(ii,ii,0,im) = nmonkap(im)
      enddo
      endif
      enddo

      end

      subroutine kap_s

      use mvariables
      use mparameters
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mKaps_s
      use mKaps

      implicit none

      include 'mpif.h'
      include 'MPI.h'
                                                     ! interaction of a
                                                     ! segment at R = a
                                                     ! and z = 0, with a
                                                     ! segment at R = b
                                                     ! and Z = c
      real*8 R,theta,z
      real*8 thetamin, thetamax, rmax, rmin, zmax, zmin, r0
      real*8 cutoff

      real*8 rn
      integer i, ii
      real*8 radio ! total radius
      real*8 x1,x2,y1, y2, z1, z2, vect
      integer Rj, Zj, j
      character*40 filename

      real*8 rounding

      integer im
      real*8 vsphere

      integer iR,iZ,itheta

      real*8 sumtemp

      real*8 h

      h = delta/2.0

      vsphere = 4.0/3.0*pi*(Kapd*delta/2.0)**3 
     & - 4.0/3.0*pi*(Kapd*delta/2.0-h)**3

      if (readkap.ne.1) then ! do not read from file

         radio = delta*dimR
         if(rank.eq.0)
     & print*,'Sphere surface calculation, readkap =', readkap

         write(filename,'(A6, I4.4,A1,I2.2, A4)')
     & 'KapsS-', dimR,'-',Kapd, '.kai'

         open(unit=112, file=filename)

         KapsV_s = 0.0               ! vector Kaps
         seed = readseed

         cutoff = Kapd*delta/2.0

         rounding = 1e-6 ! used to prevent errors due to rounding

         do ii = 1, dimR ! center of the particle    

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

               if((Zj.ge.((Kapd-1)/2+1)).or.(Zj.le.(-(Kapd-1)/2-1)))then
                print*, z, z/delta,i,ii
                stop
               endif

               if((vect.gt.(cutoff)).or.(vect.lt.(cutoff-h)))cycle ! outside sphere or inside the core

               if(Rj.gt.dimR)Rj=dimR+1 ! goes to bulk

               KapsV_s(ii, Rj, Zj) = KapsV_s(ii, Rj, Zj) + R ! incluye el jacobiano R(segmento)

            enddo ! iR
            enddo ! iZ
            enddo ! itheta

            sumtemp = 0.0
            do Rj = 1, dimR+1
            do Zj = -((Kapd-1)/2),((Kapd-1)/2)  ! Maximum z = (Kapd-1)/2
            sumtemp = sumtemp + KapsV_s(ii,Rj,Zj)
            enddo 
            enddo

! Normalize to Vsphere = 1

            do Rj = 1, dimR+1
            do Zj = -((Kapd-1)/2),((Kapd-1)/2)  ! Maximum z = (Kapd-1)/2
            KapsV_s(ii, Rj, Zj) = KapsV_s(ii, Rj, Zj)/sumtemp
            write(112,*)ii,Rj,Zj,KapsV_s(ii,Rj,Zj)
            enddo
            enddo

         end do                 ! ii -> Particle position
         close(112)
      else                      ! 
         write(filename,'(A12, I4.4,A1,I2.2, A4)')
     & '../../KapsS-', dimR,'-',Kapd, '.kai'
         open(unit=112, file=filename)
         if(rank.eq.0)print*,'Read Spheres'
         do j = 1, dimR*(dimR+1)*Kapd
            read(112,*)ii,Rj,Zj,KapsV_s(ii,Rj,Zj)

         enddo
      endif
 
!      do im = 1, N_monomer
!      Kaps(:,:,:,im) = KapsV(:,:,:)*nmonkap(im)/vsphere
!      enddo

      end


