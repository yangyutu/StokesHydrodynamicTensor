# depends on FORTRAN compilers, command to invloke the fortran compiler
FC = gfortran

#please check your compiler for these options
FFLAGS = -O3 


2d  =    deltat.o farfldres.o calwallgnsub.o wallgreenfun.o\
	      grndrm.o invert.o lubksb.o ludcmp.o \
	     main_dss.o ppexct.o ppinf.o pwexct.o pwff.o sortds.o \
           $(libs)

.f.o:;  $(FC) $(FFLAGS) -c $<

2dmc: $(2d)
	  $(FC) $(FFLAGS) -o $@ $(2d)


clean:   
	rm -f core* 2dmc *.o *.il

