# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
#   FC = ifort                    # default compiler on this machine
   FC = f90                    # default compiler on this machine
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
OBJECTS = levelAS.o masses.o alfAS.o scecor.o schrqAS.o scattlenDER.o prepot.o potgen.o AF3X3LEV.o widthAS.o 

level: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o nASlev.x

# option to clean everything up afterwards
#
#clean:
#         rm -f *.o lev.x
#--------------------------------------------------------------
# Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
