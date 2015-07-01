         subroutine chains_definitions

         use mparameters
         use mparameters_chain
         implicit none
         include 'mpif.h'
         include 'MPI.h'
         real*8 z_center
         integer i

         if(rank.eq.0) then
         print*, 'Using chain definitions for the NPC'
         print*, '20 different chains in the system'
         endif

         N_chains = 20
         ALLOCATE (long(N_chains))
         ALLOCATE (chainsperdelta(N_chains))
         ALLOCATE (zposition(N_chains))
         ALLOCATE (segtype(N_chains, maxlong)) 

         z_center   = delta*dimZ/2.0  ! NPC center z[nm]

         i=0
         chainsperdelta(:) = 8.0

         i=i+1
         ! Nup42 #1
         long(i)=382
         call read_seq(i) 
         zposition(i)      = 13 + z_center   ! in nm           

         ! Nup159 #2
         i=i+1
         long(i)=685
         call read_seq(i) 
         zposition(i)      = 13 + z_center     

         ! Nup116 #3
         i=i+1
         long(i)=789
         call read_seq(i) 
         zposition(i)      = 13 + z_center

         ! Nsp1 #4
         i=i+1
         long(i)=617
         call read_seq(i) 
         zposition(i)      = 11 + z_center
                           
         ! Nup100 #5
         i=i+1
         long(i)=800
         call read_seq(i) 
         zposition(i)      = 10 + z_center
                 
         ! Nup59 #6
         i=i+1
         long(i)=206
         call read_seq(i) 
         zposition(i)      =  8 + z_center
         
         ! Nup53 #7
         i=i+1
         long(i)=227
         call read_seq(i) 
         zposition(i)      =  7 + z_center
                         
         ! Nsp1 #8
         i=i+1
         long(i)=617
         call read_seq(i) 
         zposition(i)      =  7 + z_center
        
         ! Nup49 #9      
         i=i+1
         long(i)=251
         call read_seq(i) 
         zposition(i)      =  4 + z_center
 
         ! Nup57 #10     
         i=i+1
         long(i)=255
         call read_seq(i) 
         zposition(i)      =  4 + z_center 
                           
         ! Nup49 #11
         i=i+1
         long(i)=251
         call read_seq(i) 
         zposition(i)      = -4 + z_center
                           
         ! Nup57 #12
         i=i+1 
         long(i)=255
         call read_seq(i) 
         zposition(i)      = -4 + z_center
                         
         ! Nup145 #13   
         i=i+1
         long(i)=433
         call read_seq(i) 
         zposition(i)      =  -5 + z_center
                           
         ! Nup53 #14
         i=i+1 
         long(i)=227
         call read_seq(i) 
         zposition(i)      = -7 + z_center 
                           
         ! Nsp1 #15
         i=i+1
         long(i)=617
         call read_seq(i) 
         zposition(i)      = -7 + z_center
                           
         ! Nup59 #16
         i=i+1
         long(i)=206
         call read_seq(i) 
         zposition(i)      = -8 + z_center
         
         ! Nsp1 #17
         i=i+1
         long(i)=617
         call read_seq(i) 
         zposition(i)      = -11 + z_center
                           
         ! Nup1 #18
         i=i+1
         long(i)=857
         call read_seq(i) 
         zposition(i)      = -12 + z_center
                           
        ! Nup60 #19
         i=i+1
         long(i)=151
         call read_seq(i) 
         zposition(i)      = -14 + z_center

         ! Nup145N #20
         i=i+1
         long(i)=433
         call read_seq(i) 
         zposition(i)      = -14 + z_center

        if (i.ne.N_chains) then 
        if(rank.eq.0)print*, 'BAD NUMBER OF CHAINS for seq_type 13'
        stop
        end if                                  ! MK END ADDED

        end

C*************************************************************
      subroutine read_seq(i)

      use mparameters_chain

      implicit none
      include 'mpif.h'
      include 'MPI.h'

      character*100 filename 
      integer i, nseq, j

          write(filename,'(A15, I3.3, A4)') 'polymer_6color_',i,'.txt'  

      open(unit=2110+i,file=filename)
          
          !print*, filename, i, 'opened'
         
          read(2110+i, *), nseq
c          print*, 'nseq', nseq
          
          if (nseq.eq.long(i)) then 
                do j=1,long(i)
                        read(2110+i, *), segtype(i,j) 
                enddo
          else
                print*, 'Error in sequence file', 2110+i
                call MPI_FINALIZE(ierr) ! finaliza MPI
                stop
          endif

          close(2110+i) 
          
          return 
      end



