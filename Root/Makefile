# for C++ define  CC = g++
CC = g++
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++17
LDFLAGS = $(shell root-config --glibs) -lEG -lGeom -L./../src/ -lKMCDetFwd Loader_C.so
INCLUDE = -I$(ROOTSYS)/include -I/usr/include -I./../src/ -I./

src = $(wildcard *.C)
obj = $(src:.c=.o)

all: src load link digi reco

src:
	cd ../src; make; cd -

load:
	root -l -b -q load.C

link:
	ln -sf ../src/KMCDetFwdDict_rdict.pcm KMCDetFwdDict_rdict.pcm
	ln -sf ../src/libKMCDetFwd.so libKMCDetFwd.so

digi: Digitization.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 digi

reco: Reconstruction.C
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCLUDE)
	chmod 755 reco


.PHONY: clean
clean:
	rm -f digi reco *.o *_.cxx* *_cxx* *ACLiC* *.so dict* *.d *.pcm
	rm -rf *.dSYM
