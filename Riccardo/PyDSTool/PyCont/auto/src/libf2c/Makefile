# ===> Makefile
#
#   Makefile            Randy Paffenroth            6/1999
#
#   used to compile the f2c functions needed for the C version of AUTO
F2CSRC = cabs.c    d_lg10.c  i_nint.c  pow_di.c  r_lg10.c  z_exp.c\
         d_imag.c  d_sign.c  i_dnnt.c  pow_dd.c  pow_ii.c  z_abs.c   z_log.c

F2COBJ = cabs.o    d_lg10.o  i_nint.o  pow_di.o  r_lg10.o  z_exp.o\
         d_imag.o  d_sign.o  i_dnnt.o  pow_dd.o  pow_ii.o  z_abs.o   z_log.o

F2CSRC = $(wildcard *.c)

SOPT = -fno-common -fPIC
OPT = 
CFLAGS = $(OPT) -I../include

all: libf2c

libf2c:  $(F2COBJ)

clean:
	rm -f *.o *.a *~
