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

         N_chains = CdimZ
         ALLOCATE (long(N_chains))
         ALLOCATE (chainsperdelta(N_chains))
         ALLOCATE (zposition(N_chains))
         ALLOCATE (segtype(N_chains, maxlong)) 

         z_center =  delta*dimZ/2.0  ! NPC center z[nm]

         chainsperdelta(:) = sigma*2.0*pi*delta*CdimR*delta ! number of chains per cylindrical segment of lenght delta

         do i = 1, N_chains
         ii = i-int(N_chains/2) 

         long(i) = maxlong
         if(mod(N_chains,2).eq.1) then
             zposition(i) = z_center - ii*delta  + delta
         endif
        
         if(mod(N_chains,2).eq.0) then
            zposition(i) = z_center - ii*delta  + delta/2.0
         endif
         segtype(i, :) = 1
         segtype(i, maxlong) = 2
         enddo  
        end





