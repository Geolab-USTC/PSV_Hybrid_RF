########################################################################
#                                                                      #
# PSV Hybrid code makefile script                                      #
#                                                                      #
# Written by: Lianxing Wen, Seismological Lab.  Caltech.               #
# Last modified: Sat. 11, Aug. 1996.                                   #
#                                                                      #
########################################################################

CC   = gcc
FC   = gfortran
LIBS = -lm -L../getpar/lib -lget
GETPARINC  = -I../getpar/include
BIN  = ../bin

all: psvfd_rf aserpsvfd_rf demult clean

takeGreens: takeGreens.o readmodel.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

genmodel: genmodel.o readmodel.o interpl.o nrutil.o messy_mod.o dsyssolve.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

psvfd_rf: psvfd_rf.o readmodel.o interpl.o nrutil.o messy_mod.o dsyssolve.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

aserpsvfd_rf: aserpsvfd_rf.o aser_psv.o readmodel.o aser.rout_noplot.o dcsource.o nrutil.o
	$(FC) -o $(BIN)/$@ $^ $(LIBS)

demult: demult.o readmodel.o
	$(CC) $^ -o $(BIN)/$@ $(LIBS)

clean:
	rm -f *.o
