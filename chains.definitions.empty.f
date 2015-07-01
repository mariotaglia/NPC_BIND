         subroutine chains_definitions

         use mparameters
         use mparameters_chain
         use mvariables

         implicit none
         include 'mpif.h'
         include 'MPI.h'
         real*8 z_center
         integer i, ii

         if(rank.eq.0) then
         print*, 'Simple model for cylindrical nanopore at constant
     &  sigma'
         endif

         N_chains = 1
         ALLOCATE (long(N_chains))
         ALLOCATE (chainsperdelta(N_chains))
         ALLOCATE (zposition(N_chains))
         ALLOCATE (segtype(N_chains, maxlong)) 

         z_center =  delta*dimZ/2.0  ! NPC center z[nm]

         chainsperdelta(1) = 0.0 ! number of chains per cylindrical segment of lenght delta

         long(1) = maxlong
         zposition(1) = z_center 
         segtype(1, :) = 1
        end




