###########################################################
# fragment for ROOT compilation

ROOTCONFIG    = $(ROOTSYS)/bin/root-config
CINT          = $(ROOTSYS)/bin/rootcint

ROOTLIBS      = $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS     = $(shell $(ROOTCONFIG) --glibs) -lSpectrum
ROOTCFLAGS    = $(shell $(ROOTCONFIG) --cflags) 

CXX           = g++
LD            = g++
FOR           = 
SOFLAGS       = -shared
GLIBS         = -L/usr/X11R6/lib -lXpm -lX11 $(ROOTGLIBS)
CXXFLAGS      =  -Wall $(ROOTCFLAGS) -Iinclude/ -g
LIBS          = $(ROOTLIBS) -lMinuit

all     :	Cosmic.exe Compress.exe Pedestal.exe MIPsv2.exe DyRel.exe   
	
clean:
	rm -rf obj/* dict/* *.so;

obj/%.o : src/%.C include/Constant.h
	mkdir -p obj/; $(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.co: %.dict.cc
	mkdir -p obj/; $(CXX) $(CXXFLAGS) -c -o $@ $<; mv $< $*.dict.h dict/.

%.dict.cc: include/%.h 
	mkdir -p dict/; $(CINT) -f $@ -c -I$(ROOTSYS)/include $< 


all: obj/%.o  
	$(LD)    $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS) obj/%.o  -o $@
#electron.exe:  obj/ElecIDAna.o  obj/EMuHist.o obj/mainElectron.o obj/W_MT.o 
#	$(LD)    $(ROOTCFLAGS) $(LIBS) $(ROOTLIBS) obj/ElecIDAna.o  obj/EMuHist.o obj/mainElectron.o  obj/W_MT.o -o $@
