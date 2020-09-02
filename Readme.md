Strong field LUXE events (BPPP or Trident) propagated through a dipole magnet and hitting an ALPIDE-based tracker.
Tracks are reconstructed individually (so no combinatorics solving).
The C++ (ROOT-based) detector simulation part of this project is based on a version of ALICE software written for the NA60 experiment.
This is improved and converted for LUXE's needs here.
More details about the code are in: https://cernbox.cern.ch/index.php/s/flL3gcDe9W6pJ3V.
The reconstructed tracks are written for further analysis using python+ROOT. 


- Basic setup
  - have ROOT6 and python installed
  - setup environment for ROOT6 and python
      - `source setupROOT6.binaries.python2.sh` in my case
  - make the data and output dirs:
      - `mkdir -p data/root`
      - `mkdir -p data/stdhep/bppp`
      - `mkdir -p data/stdhep/trident`
      - `mkdir -p output/root`
      - `mkdir -p output/pdf`
  - put Tony's signal stdhep files in data/stdhep/bppp or trident dirs

- Compilation
  - hit `make`

- Convert the "truth signal" STDHEP files into ROOT files
  - open a different shell than the previous one
  - go to the python dir
  - setup python and ROOT
      - `source setupROOT6.brew.python3.sh` in my case
  - to convert the stdhep files to a ROOT TTree, run
      - `python stdhep2root.py -p bppp`
      - `python stdhep2root.py -p trident`
  - you will see 2 files in the data/root/ directory

- Generate some toy "truth background" events directly in ROOT files
  - go to the python dir
  - setup python and ROOT
      - `source setupROOT6.brew.python3.sh` in my case
  - run the background toy generation
      - `python BackgroundGeneratorTruth -p bppp -nevents 1000 -nbckgrtrk 15000`
      - `python BackgroundGeneratorTruth -p trident -nevents 100 -nbckgrtrk 20000`

- Run the "digitisation" step on the "truth signal+background" ROOT files
  - this will transform the truth tracks to a collection of clusters according to the realistic detector response
  - note that in this step (unlike in the next step) the code runs on a single track at a time!
  - look at the setup/setupLUXE_bppp.txt or setup/setupLUXE_trident.txt to understand and change the geometry and materials
  - cd Root/
  - run the signal digitisation step:
      - `root -b -q load.C 'Digitization.C+("bppp")'`
      - `root -b -q load.C 'Digitization.C+("trident")'`
  - this will produce a number of files in the data/root/ directory

- Run reconstruction step starting from the digitised clusters output:
  - cd Root/
  - run the recon step:
      - `root -b -q load.C 'runLUXEeeRecoFromClusters.C+("bppp")'`
      - `root -b -q load.C 'runLUXEeeRecoFromClusters.C+("trident")'`

- Run a basic analysis on the reconstruction step output
  - go to the python dir
  - setup python and ROOT
      - `source setupROOT6.brew.python3.sh` in my case
  - run the reco analysis:
      - `python analysis_from_clusters_reco.py -p bppp`
      - `python analysis_from_clusters_reco.py -p trident`
  - check the data and output dirs
  - run some higher-level analyses:
      - `python Resolutions.py -p bppp`
      - `python Resolutions.py -p trident`
      - `python SelectionPlots.py -p bppp`
      - `python SelectionPlots.py -p trident`
  
- Generate and visualise an event display from the reco output:
  - generate the display:
     - `python EventDisplayGenerate.py -p bppp -i 0`
     - `python EventDisplayGenerate.py -p trident -i 0`
  - visualise the display:
     - `python EventDisplayVisualize.py -p bppp`
     - `python EventDisplayVisualize.py -p trident`
  

