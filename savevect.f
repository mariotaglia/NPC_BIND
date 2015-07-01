      subroutine savevect(avvect, title, counter, counter2)



      use mlookup
      use mparameters
      use mncells

      implicit none
      real*4 varray(dimR,dimZ,2)
      real*8 avvect(dimR*dimZ,2)
      integer ix, iy, iR, i, jx, jy, jz, iZ, iC

      integer counter, counter2
      character*5 title
      character*6 titlez
      character*30 filename, tempc
      real*8 posx, posy, posz
      integer maxT, iT
      real*4 singlepres

      real*4 vectX, vectY


      varray = 0.0 ! hay que extraer las celdas con valor -1000

      do iC = 1, ncells
      varray(indexa(iC,1),indexa(iC,2),:)=avvect(iC,:)
      enddo

!!! Guarda en formato vtk


      maxT = 1
      write(filename,'(A5, A1,I3.3, A1, I3.3, A4)')title,'.',
     &  counter, '.',  counter2, '.vtk'
      open(unit=45, file=filename)
      write(45,'(A)')'# vtk DataFile Version 2.0'
      write(45,'(A)')'vectors'
      write(45,'(A)')'ASCII'
      write(45,'(A)')'DATASET STRUCTURED_GRID '
      write(45,'(A, I5, A1, I5, A1, I5)')
     & 'DIMENSIONS', dimZ+1, ' ', maxT+1, ' ',dimR+1
      write(45,'(A, I8, A)')'POINTS ',(dimz+1)*(maxT+1)
     % *(dimR+1),' float'

      do iR = 0, dimR
        do iT = 0, maxT
          do iz = 0, dimZ

      posx = sin(dfloat(iT)/maxT*2.0*3.14159)*(dfloat(iR)-0.5)*delta
      posy = cos(dfloat(iT)/maxT*2.0*3.14159)*(dfloat(iR)-0.5)*delta
      posz = iz*delta - delta/2.0

            write(45, *)
     & posx,'   ', 
     & posy
     &, '   ', posz 
           enddo
         enddo
      enddo


      write(45,'(A, I8)')'CELL_DATA ', dimR*dimZ*maxT
      tempc = 'VECTORS vectors float'
      write(45,'(A)')tempc

       do iR = 1, dimR
        do iT = 1, maxT
          do iZ = 1, dimZ

          vectX = sin(dfloat(iT)/maxT*2.0*3.14159)*varray(iR,iZ,1) 
          vectY = cos(dfloat(iT)/maxT*2.0*3.14159)*varray(iR,iZ,1)

            write(45,*)vectX, vectY, varray(iR,iZ,2)
          enddo
        enddo
       enddo
       close(45)
       end

 
