#
# COLA_DAM
#   

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
#OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc

CC      = mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP)
CFLAGS += -D_DAM_SAVEMEM #-D_LONGIDS #-D_LIGHTCONE #-D_DAM_NOMEASURE -D_DAM_NOT_WRITE_INIT
LIBS    := -lm

# Define paths of FFTW3 & GSL libraries if necessary.

FFTW3_DIR ?= /home/damonge
GSL_DIR   ?= 

DIR_PATH = $(FFTW3_DIR) $(GSL_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC = COLA_DAM
all: $(EXEC)

OBJS := src/main.o
OBJS += src/read_param.o src/msg.o src/cosmo.o
OBJS += src/pm.o src/cola.o src/comm.o src/move.o
OBJS += src/snap_io.o src/timer.o src/mem.o

LIBS += -ldl
LIBS += -lgsl -lgslcblas
LIBS += -lfftw3f_mpi -lfftw3f


ifdef OPENMP
  LIBS += -lfftw3f_omp
endif

COLA_DAM: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

src/main.o: src/main.c src/common.h
src/cola.o: src/cola.c src/common.h
src/comm.o: src/comm.c src/common.h
src/mem.o: src/mem.c src/common.h
src/move.o: src/move.c src/common.h
src/msg.o: src/msg.c src/common.h
src/pm.o: src/pm.c src/common.h
src/cosmo.o: src/cosmo.c src/common.h
src/read_param.o: src/read_param.c src/common.h
src/timer.o: src/timer.c src/common.h
src/snap_io.o: src/snap_io.c src/common.h

.PHONY: clean run dependence
clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2)

cleaner :
	rm -f $(EXEC) $(OBJS) $(OBJS2) *~ src/*~

run:
	mpirun -n 2 ./COLA_DAM param.ini

dependence:
	gcc -MM -MG *.c
