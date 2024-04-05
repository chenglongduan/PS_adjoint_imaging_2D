
include_dir1 = ${RSFROOT}/include
include_dir2 = /opt/ohpc/pub/libs/gnu7/openblas/0.2.20/include
include_dir3 = /opt/ohpc/pub/libs/gnu7/mvapich2/fftw3/3.3.8/include
link_dir1 = ${RSFROOT}/lib
link_dir2 = /opt/ohpc/pub/libs/gnu7/openblas/0.2.20/lib
link_dir3 = /opt/ohpc/pub/libs/gnu7/mvapich2/fftw3/3.3.8/lib


CC = mpiicc
CFLAGS = -Wall -traceback -O3 -x c -std=gnu99 -pedantic -fopenmp -I$(include_dir1) -I$(include_dir2) -I$(include_dir3)
LFLAGS = -L$(link_dir1) -L$(link_dir2) -lm -lmpi -lrsf -lgomp -L$(link_dir3)

# Compile source codes all at once w/o listing the names
#SRC_CODE = $(wildcard *.c)
#OBJ = $(patsubst %.c,%.o,$(wildcard *.c))

# Compile select source codes and list the names
.c.o:
	$(CC) $(CFLAGS) -c $<
SRC_CODE = bndr_rw.c \
           coef.c \
           CPML_profile.c \
           cross_correlation.c \
           main.c \
           matrix_transpose.c \
           outsnap.c \
           pad_bound.c \
           receiver.c \
           source.c \
           step_forward.c \
           outsnapimg.c
OBJ = $(SRC_CODE:%.c=%.o)

# Link
EXE = x_acrtm_ps
$(EXE) : $(OBJ)
	$(CC) -o $@ -pthread -fopenmp $(OBJ) $(LFLAGS)

clean :
	rm -rf *.o *.exe $(EXE)

all : clean $(EXE)
