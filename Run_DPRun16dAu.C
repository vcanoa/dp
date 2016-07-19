void Run_DPRun16dAu(const char *outFile = "dummy.root")
{
    gSystem->Setenv("ODBCINI","/opt/phenix/etc/odbc.ini.test");	    
    gSystem->Load("libDPRun16dAu.so");
    gSystem->Load("libTOAD");
    
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);
    
    TOAD toad_loader("PhotonConversionAnalysis");
    toad_loader.SetVerbosity(4);
    string lookupfile_location = toad_loader.location("lookup_3D.root");
    cout << "Loading TOAD" << endl;
    
    mDstToPhotonEvent *rd = new mDstToPhotonEvent(outFile, lookupfile_location.c_str());
    se->registerSubsystem( rd );
    cout << "Module attached" << endl;
}
void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  indata.push_back("DST_EVE");
  return;
}
