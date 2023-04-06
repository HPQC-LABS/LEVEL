#   as in:   make FC=ifort
ifndef $FC
#  FC = nagfor
#  FC = f90                
   FC = gfortran          
endif
#
#   ad in:   make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
  else                           
    ifeq ($(debug),2)
        FFLAGS =  -C -g       # debug=2  (higher-level debug option)
    endif
  endif
OBJECTS = levelAS.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o level.x

# To run the code, execute:

