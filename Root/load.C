{
  cout << "Soft-linking the KMCDetFwd library to here..." << endl;
  gSystem->Exec("ln -sf ../src/libKMCDetFwd.so libKMCDetFwd.so");
  cout << "Loading the KMCDetFwd library..." << endl;
  gSystem->Load("libKMCDetFwd.so");
  cout << "Done..." << endl;
}
