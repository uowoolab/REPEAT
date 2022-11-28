#
#  Makefile for Sies2xcf set of tools
#
# Compiler and options:
#  FC     = ifort
#  FFLAGS = -mp1 -O3 -pc80 -prec_div -r8 -w
#FC=ifort -O2
#FFLAGS=-g -fixedform -static -m64 \
	   -TENV\:simd_zmask=OFF -TENV\:simd_imask=OFF \
	   -TENV\:simd_omask=OFF -O3 -OPT\:Ofast \
	   -fno-math-errno -I.
FC=ifort
FFLAGS= -g
LFLAGS = $(FFLAGS)
#
EXECS = eig2bxsf xv2xsf md2axsf rho2xsf vib2xsf siesta2repeat
.f.o      : $(FC) -c $(FFLAGS) $<

all: $(EXECS)
   eig_obj = eig2bxsf.o inver3.o  opnout.o  
   xv_obj  = xv2xsf.o opnout.o 
   md_obj  = md2axsf.o test_md.o makebox.o fillbox.o inver3.o  hit.o \
             wraxsf1.o wraxsf2.o opnout.o test_ani.o
   rho_obj = rho2xsf.o read_xv.o makebox.o fillbox.o inver3.o  hit.o \
             intpl04.o opnout.o   
   siesta2repeat_obj = siesta2repeat.o read_xv.o makebox.o fillbox.o inver3.o  hit.o \
             intpl04.o opnout.o   
   vib_obj = vib2xsf.o read_xv.o makebox.o fillbox.o inver3.o  hit.o \
	     displa.o  read_ev.o itochar.o w_arrow.o w_movie.o opnout.o 

 eig2bxsf : $(eig_obj)
	$(FC) $(LFLAGS) -o $@ $(eig_obj)
  xv2xsf   : $(xv_obj)
	$(FC) $(LFLAGS) -o $@ $(xv_obj)
  md2axsf  : $(md_obj)
	$(FC) $(LFLAGS) -o $@ $(md_obj)
  rho2xsf  : $(rho_obj)
	$(FC) $(LFLAGS) -o $@ $(rho_obj)
  siesta2repeat  : $(siesta2repeat_obj)
	$(FC) $(LFLAGS) -o $@ $(siesta2repeat_obj)
  vib2xsf  : $(vib_obj)
	$(FC) $(LFLAGS) -o $@ $(vib_obj)

	

clean:
	@echo " Cleaning up"
	rm -f $(EXECS) *.o core

