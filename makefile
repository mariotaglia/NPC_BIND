#Unix makefile for fortran-file	

# put the name of the target program here
TARGET = nanoporo
SRC = modules.f nanoporo.f cadenas.f encode.f fkfun.f  kinsol.f  monomers.definitions.f  sphere.f chains.definitions-simple.f  factorcurv.f  free_energy.f  lookup.f rands.f  globals.f  saveresults.f kais.f savevect.f

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 #-fbounds-check ${F90FLAGS}

LDFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
CFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

 #${LINK_F90_STATIC}

FF = mpif77 #${F90}


all:	$(TARGET)

$(TARGET): $(SRC:.f=.o)
	$(FF) -o $(TARGET) $(SRC:.f=.o) $(LFLAGS) $(LDFLAGS) $(CFLAGS)
	mv nanoporo ~/bin/NPC_BIND/
.f.o:
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(LDFLAGS) $(CFLAGS)

install: all
	cp $(TARGET) ~/bin


clean:	
	@rm -f $(SRC:.f=.o) $(SRC:.f=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif


















































