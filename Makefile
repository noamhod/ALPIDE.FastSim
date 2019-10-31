# for C++ define  CC = g++
CC = g++
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++17
LFLAGS = $(shell root-config --glibs) -lEG -lGeom
INC =	-I$(ROOTSYS)/include -I/usr/include  -I./
TGT =	libKMCDetFwd.so
DICT=	KMCDetFwdDict.cxx
DICTO=	KMCDetFwdDict.o

SRC =  TTreeStream.cxx GenMUONLMR.cxx KMCClusterFwd.cxx  KMCDetectorFwd.cxx  KMCFlukaParser.cxx  KMCLayerFwd.cxx  KMCProbeFwd.cxx  KMCUtils.cxx NaMaterial.cxx AliLog.cxx TrackPar.cxx

HDR =	$(SRC:.cxx=.h) 

OBJ = 	$(SRC:.cxx=.o)


.PHONY: depend clean

all: 	$(TGT)
	@echo creating libKMCDetFwd.so

$(TGT):	$(OBJ) $(DICTO)
	$(CC) $(CFLAGS)  -shared -o $(TGT) $(OBJ) $(DICTO) $(shell root-config --ldflags) $(LFLAGS)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

%.o : %.cxx
	$(CC) $(CFLAGS) $(INC) -c $<  -o $@

clean:
	rm -rfv $(DICT) *.o *~ *.so *.pcm *.d *_.cxx

$(DICT): $(HDR) KMCDetFwfLinkDef.h
	rootcint -f $@  $(INC) $(HDR) $^


depend: $(SRC)
	makedepend $(INC) $^

# DO NOT DELETE THIS LINE -- make depend needs it
