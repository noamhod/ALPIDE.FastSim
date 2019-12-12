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

- Run the reconstruction on the "truth" ROOT files
  - look at the setup/setupLUXE.txt to understand/change the geometry and materials
  - in runLUXEeeReco.C, find this line: `TString process = "trident";  /// trident or bppp`
  - change (hardcode) as needed between "bppp" and "trident"
  - run `root -b -q load.C runLUXEeeReco.C+`
  - repeat the 2 steps above to flip between "bppp" and "trident"
  - this will produce a number of files in the data/root/ directory

- Analyse the reconstructed output:
  - go to the python dir
  - setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  - run the analysis:
     - `python analysis.py -p bppp`
     - `python analysis.py -p trident`