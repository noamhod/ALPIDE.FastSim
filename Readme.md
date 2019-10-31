Strong field LUXE events (BPPP or Trident) propagated through a dipole magnet and hitting an ALPIDE-based tracker.
Tracks are reconstructed individually (so no combinatorics solving).
The C++ (ROOT-based) detector simulation part of this project is based on a version of ALICE software written for the NA60 experiment.
This is improved and converted for LUXE's needs here.
More details about the code are in: https://cernbox.cern.ch/index.php/s/flL3gcDe9W6pJ3V.
The reconstructed tracks are written for further analysis using python+ROOT.

- Basic setup
  [1] have ROOT6 and python installed
  [2] setup environment for ROOT6 and python (`source setupROOT6.binaries.python2.sh` in my case)
  [3] put Tony's stdhep files in data/stdhep/bppp and data/stdhep/trident
  [4] hit `make`

- Convert the "truth" STDHEP files into ROOT files
  [1] go to the python dir
  [2] setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  [3] to convert the stdhep files to a ROOT TTree, run
      - python stdhep2root.py -p bppp
      - python stdhep2root.py -p trident
  [4] you will see 2 files in the data/root/ directory

- Run the reconstruction on the "truth" ROOT files
  [1] look at the setup/setupLUXE.txt to understand/change the geometry and materials
  [2] in runLUXEeeReco.C, find this line: `TString process = "trident";  /// trident or bppp`
  [3] change (hardcode) as needed between "bppp" and "trident"
  [4] run `root -b -q load.C runLUXEeeReco.C+`
  [5] repeat (2) and (3) for "bppp" or "trident"
  [6] this will produce a number of files in the data/root/ directory

- Analyse the reconstructed output:
  [1] go to the python dir
  [2] setup python and ROOT (`source setupROOT6.brew.python3.sh` in my case)
  [3] run the analysis:
     - `python analysis.py -p bppp`
     - `python analysis.py -p trident`