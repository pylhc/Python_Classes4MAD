all: transport.so

%.so: %.f90 Makefile
	f2py -c --opt=-O3 --f90flags=-Wall --fcompiler=gnu95 $< -m $*

clean:
	-rm *.so
