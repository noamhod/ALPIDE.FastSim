# for C++ define  CC = g++
CC = g++
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++17
LDFLAGS = $(shell root-config --glibs) -lEG -lGeom -L./../src/ -lKMCDetFwd Loader_C.so
### depends on where the root header file is located
#INCLUDE = -I$(ROOTSYS)/include -I/usr/include -I./../src/ -I./
INCLUDE = -I$(ROOTSYS)/include -I$(ROOTSYS)/include/root -I/usr/include -I./../src/ -I./

src = $(wildcard *.C)
obj = $(src:.c=.o)

all: src load link digi digiE320 digiNew reco recoNew recoFull

src:
	cd ../src; make; cd -

load:
	root -l -b -q load.C

link:
	ln -sf ../src/KMCDetFwdDict_rdict.pcm KMCDetFwdDict_rdict.pcm
	ln -sf ../src/libKMCDetFwd.so libKMCDetFwd.so

digiE320: DigitizationE320.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 digiE320

digi: Digitization.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 digi

digiNew: DigitizationNew.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 digiNew

reco: Reconstruction.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 reco

recoNew: ReconstructionNew.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 recoNew

recoFull: ReconstructionFullSim.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 recoFull


.PHONY: clean
clean:
	rm -f digi digiE320 reco recoNew *.o *_.cxx* *_cxx* *ACLiC* *.so dict* *.d *.pcm
	rm -rf *.dSYM
