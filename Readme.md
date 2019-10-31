
Builing the library requires ROOT installed and sources (Makefile assumes
ROOT 6, but 5 also should work but some compiliation options may require change)

Build the library by calling ``make``

One can run a simple test: muons from omega->mu+ mu- decays at {0,0,0} propagated
through the dipole field towards 5 silicon planes in field-free region.

Run the test as:

``root -b -q load.C runDimuTest.C+``

The output for every generated omega is streamed into the root tree, storing
the generated and reconstructed single track and dimuons kinematics
(as a TLorentzVector), number of muon traks reconstructed as well as their
fit chi2's etc.
````
root [1] auto* tr = (TTree*)gFile->Get("res");
root [2] tr->Draw("rec.M() - gen.M()","nrec==2"); // difference in rec. and gen. masses
````
Very incomplete description of the code: https://cernbox.cern.ch/index.php/s/flL3gcDe9W6pJ3V


Basic setup
1. have ROOT6 and python installed
2. setup environment for ROOT6 and python (`source setupROOT6.binaries.python2.sh` in my case)
3. put Tony's stdhep files in data/stdhep/bppp and data/stdhep/trident
4. hit `make`

Convert the "truth" STDHEP files into ROOT files
1. go to the python dir
2. setup python3 and ROOT6 (`source setupROOT6.brew.python3.sh` in my case)
3. to convert the stdhep files to a ROOT TTree, run
   - python stdhep2root.py -p bppp
   - python stdhep2root.py -p trident
4. you will see 2 files in the data/root/ directory

Run the reconstruction on the "truth" ROOT files
1. in runLUXEeeReco.C, find this line: TString process = "trident";  /// trident or bppp
2. change (hardcode) as needed between "bppp" and "trident"
3. run `root -b -q load.C runLUXEeeReco.C+`
2. repeat (2) and (3) for "bppp" or "trident"
4. this will produce a number of files in the data/root/ directory

Analyse the 
1. go to the python dir
2. setup python3 and ROOT6 (`source setupROOT6.brew.python3.sh` in my case)
3. run 