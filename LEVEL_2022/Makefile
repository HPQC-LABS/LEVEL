#   as in:   make FC=ifort
ifndef $FC
#  FC = nagfor
#  FC = f90                
   FC = gfortran          
endif

#   as in:   make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
  else                           
    ifeq ($(debug),2)
        FFLAGS =  -C -g       # debug=2  (higher-level debug option)
    endif
  endif

OBJECTS = levelAS.o
level: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o level.x

# To run the code, execute:
# ./level.x < input.5 > fort.6

# [nike@gra-login3 LEVEL_2022]$ make
# ifort  -O2 -u              -c -o levelAS.o levelAS.f
# ifort  -O2 -u             levelAS.o  -o level.x
# 
# [nike@gra-login1 LEVEL_2022]$ make DEBUG=1
# ifort  -g -u               -c -o levelAS.o levelAS.f
# ifort  -g -u              levelAS.o  -o level.x
# 
# nike@gra-login1 LEVEL_2022]$ make DEBUG=2
# ifort  -C -g -u         -c -o levelAS.o levelAS.f
# ifort  -C -g -u        levelAS.o  -o level.x
# 
# Similar to DEBUG=1 but with -O0 (could try -O1 and -O2 too):
# ifort -g -O0 -u -c -o levelAS.o levelAS.f
# ifort -g -O0 -u  levelAS.o -o level.x
