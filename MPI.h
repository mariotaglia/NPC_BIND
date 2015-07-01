! MPI

       integer rank, size, ierr

       integer stat(MPI_STATUS_SIZE)
       integer tag, source, dest
       parameter(tag = 0)
       integer err
       integer flagsolver

       integer ier_tosend
       double  precision norma_tosend

       common /MPI/ rank, size, ierr, flagsolver

