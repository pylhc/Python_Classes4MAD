# ===> Makefile
#
#   Makefile            Drew LaMarh            4/2006
#
#   used to compile the module functions needed for the PyDSTool version of AUTO
#
SOPT = -fno-common -fPIC
OPT = -Wno-long-double -DPYTHON_MODULE

AUTOMODSRC = automod.c interface.c

AUTOMODOBJ = automod.o interface.o

CFLAGS = $(OPT) -I../src/include -Iinclude -I/sw/include/python2.3

LIBDIR = -L/sw/lib -L../lib
LDFLAGS = -dynamiclib -undefined dynamic_lookup -install_name $(AUTO_DIR)/lib/libautomod.dylib
LIBS = 

all: swig libautomod

libautomod:  $(AUTOMODOBJ)
	#gcc $(LIBDIR) $(LDFLAGS) $(LIBS) $(AUTOMODOBJ) -o libautomod.dylib
	#mv libautomod.dylib ../lib
	ar rc libautomod.a $(AUTOMODOBJ)
	mv libautomod.a ../lib

swig:	
	swig -python -o automod_wrap.c automod.i

clean:
	rm -f *.o *.a *.so *~ *.pyc automod_wrap.c
	cd include;rm -f *.h~
