# This is the makefile for diablo.
# To compile the code, just type make (if environmental variables are set
# appropriately) -- otherwise you must run the go scripts.  
# Such an approach makes recompilation of the code easy, recompiling 
# only as necessary to account for recent changes to the code.
#

# As the user, set the following definitions
# Fortran 90 complier to use:
COMPILER = ifort 
#COMPILER = xlf
#COMPILER = mpif90

# Any desired options for the compiler (e.g. -O2, -g, etc.)
USEROPTS = -O2
#USEROPTS = -O3
#USEROPTS = -O2 -warn all  -fpe0 -traceback -debug inline_debug_info

cgl: cgl.f  var.o grid_def 
	$(COMPILER) $(USEROPTS) cgl.f -o cgl var.o
	

var.o:  var.f grid_def
	$(COMPILER) $(USEROPTS) -c var.f

clean:
	rm -f cgl *.o fort.* *~ *.mod *.tec ;
#	clear;

