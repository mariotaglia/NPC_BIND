      subroutine monomer_definitions

      use mparameters_monomer
      use mKaps

      implicit none

      include 'mpif.h'
      include 'MPI.h'

      integer im
      integer ios, j

      N_poorsol = 1 ! number of different kais
      N_monomer = 4

      ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
      ALLOCATE (pKbind(N_monomer), Kbind(N_monomer), Kbind0(N_monomer)
     & , rbind(N_monomer))
      ALLOCATE (henergy(N_poorsol))
      ALLOCATE (center(N_monomer))

      ios = 0
      read(nmonkapbuffer, *, iostat=ios) (nmonkap(j),j=1,N_monomer)

      do im = 1, N_monomer
      if(rank.eq.0)print*, 'Monomers type', im,' in Kap:', nmonkap(im)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      st_matrix(1,1)=1.0

      henergy(1) =  st_matrix(1,1)

! Segment type 1: neutral, hydrophilic, not bind

      zpol(1) = 0
      hydroph(1) = 0
      pKa(1) = -1.0

      rbind(1) = 0
      pKbind(1) = -1.0 

      center(1) = 0

! Segment type 2: neutral, hydrophilic, bind

      zpol(2) = 0
      hydroph(2) = 0
      pKa(2) = -1.0

      rbind(2) = 1
      pKbind(2) = pKbindread

      center(2) = 0

! type 3: positive, hydrophilic, not bind

      zpol(3) = 1
      hydroph(3) = 0
      pKa(3) = 14
      
      rbind(3) = 0
      pKbind(3) = -1.0

      center(3) = 0
! type 4: center, neutral, hydrophobic, not bind

      zpol(4) = 1
      hydroph(4) = 1
      pKa(4) = 14

      rbind(4) = 0
      pKbind(4) = -1.0

      center(4) = 1



      end

