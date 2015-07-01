         subroutine decode_one(displ, displ_one)
c Decode a 8 bit integer (displ) to a 8 size binary array (displ_one)
         implicit none
         integer*1 displ, k
         integer*1 displ_one(0:7)
         do k = 0, 7
         displ_one(k) = iand(ishft(displ,k-7),1)
         enddo
         end

         subroutine encode_one(displ, displ_one)
c Encode a 8 bit integer (displ) from a 8 size binary array (displ_one)
         implicit none
         integer*1 displ, k
         integer*1 displ_one(0:7), temp(0:7)
         displ = 0
         do k = 0,7
         temp = ishft(displ_one,7-k)
         displ = displ+temp(k)
         enddo
         end

         subroutine encode(displacement, binary, long1)
c Encodes a displacement matrix of length "long, 2" into a binary matrix of
c length int(long/2)


         use mparameters
         use mparameters_chain
         use mparameters_monomer

         implicit none
         include 'mpif.h' ! MPI libraries
         include 'MPI.h' ! MPI libraries


         integer*1 displ
         integer*1 displ_one(0:7)
         integer j, i
         integer*1 displ_temp         
         integer long1

         integer*1 displacement(maxlong+1,2)
         integer*1 binary(int(maxlong/2))

          do j = 2, long1

          select case (displacement(j, 1)) ! displacement in R
          case (0) ! no displacement
          displ_one(1+4*mod(j,2)) = 0
          displ_one(2+4*mod(j,2)) = 0
          case (1) ! shift right
          displ_one(1+4*mod(j,2)) = 0
          displ_one(2+4*mod(j,2)) = 1
          case (-1) ! shift left
          displ_one(1+4*mod(j,2)) = 1
          displ_one(2+4*mod(j,2)) = 0
          case default
          print*, 'Error in creador R, lseg should be < delta'
          print*, j, displacement(j, 1)
          call MPI_FINALIZE(ierr) ! finaliza MPI
          stop

          endselect

          select case (displacement(j, 2)) ! displacement in Z

          case (0) ! no displacement
          displ_one(0+4*mod(j,2)) = 0
          displ_one(3+4*mod(j,2)) = 0
          case (1) ! shift bottom
          displ_one(0+4*mod(j,2)) = 0
          displ_one(3+4*mod(j,2)) = 1
          case (-1) ! shift up
          displ_one(0+4*mod(j,2)) = 1
          displ_one(3+4*mod(j,2)) = 0
          case default
          print*, 'Error in creador Z, lseg should be < delta'
          call MPI_FINALIZE(ierr) ! finaliza MPI
          stop

          endselect

          call encode_one(displ_temp,displ_one)
          binary(int(j/2)) = displ_temp
          end do  ! j
        
         end

         subroutine decode(displacement, binary, long1)


         use mparameters
         use mparameters_chain
         use mparameters_monomer

         implicit none

c Decodes a displacement matrix of length "long+1, 2" from a binary matrix of
c length int(long/2)
c (uses long+1 to prevent errors in the case of even number of segments!)

         integer*1 displ
         integer*1 displ_one(0:7)
         integer j, i
         integer long1

         integer*1 displacement(maxlong+1,2)
         integer*1 binary(int(maxlong/2))

         displacement(1,1) = 0 ! R
         displacement(1,2) = 0 ! Z

         do i=1, int(long1/2)

         displ = binary(i)

         call decode_one(displ, displ_one)

         j = i * 2

         displacement(j,1) =  
     & - displ_one(1) + displ_one(2)  ! R
         displacement(j,2) = 
     & - displ_one(0) + displ_one(3)  ! Z

         j = i * 2 + 1
         
         displacement(j,1) =
     & - displ_one(5) + displ_one(6)  ! R
         displacement(j,2) =
     & - displ_one(4) + displ_one(7)  ! Z

         enddo

         end


