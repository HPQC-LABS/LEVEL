# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
FC=nagfor
#FC=gfortran
# FC=ifort 
# FC=f90    
endif
#
#option to choose level of optimization:  make DEBUG=1
ifndef DEBUG
    FFLAGS =  -O2 -u            # DEBUG not specified, so optimize
  else                           
    FFLAGS =  -g -u             # DEBUG=1  (basic DEBUG option)
    ifeq ($(DEBUG),2)           
        FFLAGS =  -C -g -u       # DEBUG=2  (higher-level DEBUG option)
    endif
#                 can add more options for  DEBUG=3 , etc., as desired
endif
#
# as usual, list the objects
#
OBJECTS = levelAS.o 

level: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o level.x

# option to clean everything up afterwards
#
#clean:
#         rm -f *.o level.x
#--------------------------------------------------------------
# Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
# To run the program do:
# ./level.x < input.5 > fort.6
#
#














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
