Strong field LUXE events (BPPP or Trident) propagated through a dipole magnet and hitting an ALPIDE-based tracker.
Tracks are reconstructed individually (so no combinatorics solving).
The C++ (ROOT-based) detector simulation part of this project is based on a version of ALICE software written for the NA60 experiment.
This is improved and converted for LUXE's needs here.
More details about the code are in: https://cernbox.cern.ch/index.php/s/flL3gcDe9W6pJ3V.
The reconstructed tracks are written for further analysis using python+ROOT. 


- Basic setup
  - have ROOT6 and python installed
  - setup environment for ROOT6 and python (`source setupROOT6.binaries.python2.sh` in my case)
  - put Tony's stdhep files in data/stdhep/bppp and data/stdhep/trident
  - hit `make`

- Convert the "truth" STDHEP files into ROOT files
  - go to the python dir
  - setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  - to convert the stdhep files to a ROOT TTree, run
      - python stdhep2root.py -p bppp
      - python stdhep2root.py -p trident
  - you will see 2 files in the data/root/ directory

- Run the "digitisation" part on the "truth" ROOT files
  - this will transform the truth tracks to a collection of clusters according to the realistic detector response
  - note that the code runs on a single track at a time!
  - look at the setup/setupLUXE_bppp.txt or setup/setupLUXE_trident.txt to understand/change the geometry and materials
  - run the signal digitisation step:
    - root -b -q load.C 'runLUXEeeReco.C+("bppp")'
    - root -b -q load.C 'runLUXEeeReco.C+("trident")'
  - this will produce a number of files in the data/root/ directory

- Basic analysis of the digitisation output:
  - in a separate shell, go to the python dir
  - setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  - run the simple analysis on the digitisation step:
    - python analysis.py -p bppp
	 - python analysis.py -p trident
  - check the data and output dirs...

- Run a "seeding demo" starting from the "digitisation" output:
  - in a separate shell, go to the python dir
  - setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  - run the demo alg and plotting script:
    - python SeedingDemo.py -p bppp
	   - python SeedingDemoPlot.py -p bppp
	 - python SeedingDemo.py -p trident
	   - python SeedingDemoPlot.py -p trident
  - check the data and output dirs...

- Run a "reconstruction demo" starting from the "seeding demo" output:
  - run the recon step:
    - root -b -q load.C 'runLUXEeeRecoFromSeeds.C+("bppp")'
	 - root -b -q load.C 'runLUXEeeRecoFromSeeds.C+("trident")'
  - in a separate shell, go to the python dir
  - setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  - run the demo plotting script:
    - python RecoDemoPlot.py -p bppp
	 - python RecoDemoPlot.py -p trident
  - check the data and output dirs...


