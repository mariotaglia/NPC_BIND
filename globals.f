      subroutine globals

      use mparameters
      use mparameters_monomer
      use mvariables
      use mprotein
      use mKaps

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries


! Input related variables
      character (len=100)  buffer,label
      integer pos
      integer, parameter :: fh = 15
      integer ios 
      integer line, linemax
      integer i, j

      character(len=50) :: filename = 'DEFINITIONS.txt'


! Control file variables

      line = 0
      ios = 0

      open(fh, file=filename)

      if(rank.eq.0)print*, 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.
     
     
      do while (ios == 0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1

! Find the first instance of whitespace.  Split label and data.
 

        pos = scan(buffer, ' ')
 
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)

       case ('pKbind')
           read(buffer, *, iostat=ios) pKbindread
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('hst')
           read(buffer, *, iostat=ios) hst
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('hst2')
           read(buffer, *, iostat=ios) hst2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('Ptype')
           read(buffer, *, iostat=ios) Ptype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('hrange')
           read(buffer, *, iostat=ios) hrange
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('saveindex')
           read(buffer, *, iostat=ios) saveindex
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('poretype')
           read(buffer, *, iostat=ios) poretype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('cadenastype')
           read(buffer, *, iostat=ios) cadenastype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dimR')
           read(buffer, *, iostat=ios) dimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dimZ')
           read(buffer, *, iostat=ios) dimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('delta')
           read(buffer, *, iostat=ios) delta
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('lseg')
           read(buffer, *, iostat=ios) lseg
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('cuantas')
           read(buffer, *, iostat=ios) cuantas
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('CdimR')
           read(buffer, *, iostat=ios) CdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('CdimZ')
           read(buffer, *, iostat=ios) CdimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('CdimZmin')
           read(buffer, *, iostat=ios) CdimZmin
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('BdimR')
           read(buffer, *, iostat=ios) BdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('TdimR')
           read(buffer, *, iostat=ios) TdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('RdimZ')
           read(buffer, *, iostat=ios) RdimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('Curvature')
           read(buffer, *, iostat=ios) curvature
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('savetodisk_type')
           read(buffer, *, iostat=ios) savetodisk_type
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('save_every')
           read(buffer, *, iostat=ios) save_every
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dseg')
           read(buffer, *, iostat=ios) dseg
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('nearbonds')
           read(buffer, *, iostat=ios) nearbonds
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('mcube')
           read(buffer, *, iostat=ios) mcube
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('ncha_max')
           read(buffer, *, iostat=ios) ncha_max
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('wantedCPUsecs')
           read(buffer, *, iostat=ios) wantedCPUsecs
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('calq')
           read(buffer, *, iostat=ios) calq
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('qprob')
           read(buffer, *, iostat=ios) qprob0
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('csalt')
           read(buffer, *, iostat=ios) csalt
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pHbulk')
           read(buffer, *, iostat=ios) pHbulk
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('infile')
           read(buffer, *, iostat=ios) infile
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('st')
         read(buffer, *, iostat=ios) st
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)


        case ('savetodisk_flag')
           read(buffer, *, iostat=ios) savetodisk_flag
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('sigma')
           read(buffer, *, iostat=ios) sigma
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('maxlong')
           read(buffer, *, iostat=ios) maxlong
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('MCsteps')
           read(buffer, *, iostat=ios) MCsteps
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('readkai')
           read(buffer, *, iostat=ios) readkai
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('readkap')
           read(buffer, *, iostat=ios) readkap
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

          case ('nxkapbulk')
           read(buffer, *, iostat=ios) nxkapbulk
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
           read(fh, *),(xkapbulks(i), i=1, nxkapbulk)
         if(rank.eq.0)
     &   print*, 'xkapbulks: ', (xkapbulks(j), j=1,nxkapbulk)

        case ('Kapd')
           read(buffer, *, iostat=ios) Kapd
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         if(mod(Kapd,2).eq.0) then
           if(rank.eq.0)print*,'Kapd should be odd!'
           stop 
         endif

         case ('nmonkap')
         nmonkapbuffer = buffer
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
          
         case ('Ppos')
          read(buffer, *, iostat=ios) Ppos
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('PdimR')
          read(buffer, *, iostat=ios) PdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('PdimH')
          read(buffer, *, iostat=ios) PdimH
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_vol')
          read(buffer, *, iostat=ios) prot_vol
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_q')
          read(buffer, *, iostat=ios) prot_q
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_q2')
          read(buffer, *, iostat=ios) prot_q2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('weakP')
          read(buffer, *, iostat=ios) weakP
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('pKaP')
          read(buffer, *, iostat=ios) pKaP
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('readseed')
          read(buffer, *, iostat=ios) readseed
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('monomerz')
          read(buffer, *, iostat=ios) monomerz
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('monomerpKa')
          read(buffer, *, iostat=ios) monomerpKa
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        end select

        endif 
        enddo     

        ncha_max=cuantas/1000+200
       if(rank.eq.0)print*,'Set ncha_max = ', ncha_max

        b2 = lseg*lseg
       if(rank.eq.0)print*,'Set b2 = ', b2
        d2 = dseg*dseg
       if(rank.eq.0)print*,'Set d2 = ', d2

       end

       subroutine array_alloc

      use mparameters
      use mchains
      use mlookup
      use mkai
      use mporesystem
      use msegme
      use mvariables
      use mkinsol
      use mparameters_chain
      use mKaps
      use mKaps_s
      use mparameters_monomer
      use posmk
      use mprotein

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      integer listkap

      listkap = 1000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ALLOCATE ARRAYS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

! chains

!        ALLOCATE (fs(cuantas))        
        ALLOCATE (inc(cuantas, maxlong))        
!        ALLOCATE (displ(cuantas, int(maxlong/2)))
        ALLOCATE (endtoend(cuantas))

! lookup

      ALLOCATE (rp(dimR*dimZ))
      ALLOCATE (rm(dimR*dimZ))
      ALLOCATE (zp(dimR*dimZ))
      ALLOCATE (zm(dimR*dimZ))
      ALLOCATE (indexa(dimR*dimZ,2)) ! 1 = R, 2 = Z

! kai
     
      ALLOCATE (Xu(dimR, dimR+1, -1:1))
      ALLOCATE (nXu(dimR*dimZ))
      ALLOCATE (Xulist_cell(dimR*dimZ, 10))
      ALLOCATE (Xulist_value(dimR*dimZ, 10))

! Kaps

      if(rank.eq.0) then

      ALLOCATE (KapsV(dimR, dimR+1, -(Kapd-1)/2:(Kapd-1)/2))
      ALLOCATE (nKap(dimR*dimZ))
      ALLOCATE (Kaplist_cell(dimR*dimZ, listkap))

      ALLOCATE (Kaps(dimR,dimR+1,-(Kapd-1)/2:(Kapd-1)/2,N_monomer))
      ALLOCATE (Kaplist_value(dimR*dimZ, listkap, N_monomer))
      ALLOCATE (kapmask(dimR*dimZ))
    
      endif

      ALLOCATE (xtotalbulk(N_poorsol))
      ALLOCATE (fdisbulk(N_monomer))
      ALLOCATE (xpotbulk(N_monomer))

! Kaps S

      if(rank.eq.0) then
      ALLOCATE (KapsV_s(dimR, dimR+1, -(Kapd-1)/2:(Kapd-1)/2))
      ALLOCATE (nKap_s(dimR*dimZ))
      ALLOCATE (Kaplist_s_cell(dimR*dimZ, listkap))
      ALLOCATE (v0r(dimR*dimZ, listkap))
      ALLOCATE (v0z(dimR*dimZ, listkap))
      ALLOCATE (Kaplist_s_value(dimR*dimZ, listkap))
      endif

! pore.system

      ALLOCATE (matriz(0:dimR+1, 0:dimZ+1))

! segme

      ALLOCATE (in1(maxlong,3))
      ALLOCATE (celda(maxlong))

! variables

      ALLOCATE (avpol(N_monomer, N_chains+2, dimR*dimZ))
      ALLOCATE (qtot(dimR*dimZ))
      ALLOCATE (receptor(dimR*dimZ))
      ALLOCATE (rhokap(dimR*dimZ))
      ALLOCATE (rhokapb(dimR*dimZ))
      ALLOCATE (psi2(dimR*dimZ))
      ALLOCATE (xsol(dimR*dimZ))
      ALLOCATE (xtotal2(N_poorsol, dimR*dimZ+2))
      ALLOCATE (xpos(dimR*dimZ))
      ALLOCATE (xneg(dimR*dimZ))
      ALLOCATE (xHplus(dimR*dimZ))
      ALLOCATE (xOHmin(dimR*dimZ))

      ALLOCATE (fdis(N_monomer, dimR*dimZ))
      ALLOCATE (funbound(N_monomer, dimR*dimZ))
      ALLOCATE (fbound(N_monomer, dimR*dimZ, listkap))
      ALLOCATE (fdisP(dimR*dimZ))
      ALLOCATE (newcuantas(N_chains))

! kinsol

      ALLOCATE (pp((2+N_poorsol)*dimR*dimZ))

! posmk

      ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
      ALLOCATE (current(maxlong, 3))
      ALLOCATE (nextbead(maxlong))

! protein

      ALLOCATE (protein(dimR, dimZ))
      ALLOCATE (proteinq(dimR, dimZ))
      ALLOCATE (proteinh(dimR, dimZ))
      ALLOCATE (proteinC(dimR*dimZ))
      ALLOCATE (proteinqC(dimR*dimZ))
      ALLOCATE (proteinhC(dimR*dimZ))

      end




