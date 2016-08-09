#ifndef __MPHOTONEVENTQA_H__
#define __MPHOTONEVENTQA_H__

#include <TFile.h>

#include "THmulf.h"
#include "SubsysReco.h"

using namespace std;

class TFile;
class TH1F;
class TH2F;
class THmulf;

class Fun4AllServer;
class CentralityReco;
class recoConsts;

class PHCompositeNode;
class PHGlobal;
class TrigLvl1;
class PHCentralTrack;
class emcClusterContainer;
class emcClusterContent;
class ReactionPlaneObject;
class ReactionPlaneSngl;
class RunHeader;

class CNTDE;
class CNTE;
class EMCC;
class PhotonEvent;


class mPhotonEventQA: public SubsysReco
{
 public:
  mPhotonEventQA(const char *outfile = "TEST.root");
  virtual ~mPhotonEventQA();
  
  int   Init(PHCompositeNode *topNode);
  int   InitRun(PHCompositeNode *topNode);
  int   process_event(PHCompositeNode *topNode);
  int   ResetEvent(PHCompositeNode *topNode) {return 0;};
  int   Reset(PHCompositeNode *topNode) {return 0;};
  int   End(PHCompositeNode *topNode);
  int   GetNodes(PHCompositeNode *topNode);
  bool  trigger_selection();
  bool  cut_for_calibration(int itrk_reco);
  void  calibrations();
  void  DC_map();
  void  Emcal_map();
  // void  RICH_map();
  void  run_by_run();
  void  Print(const char *what) const {return;};
  
 protected:   
  TFile* f;
  recoConsts* rc; 
  CentralityReco* centReco;
  Fun4AllServer* se;
  
  // Reconstruction reco;
  
  PHGlobal* phg;
  TrigLvl1* _Trig_ptr;    
  PHCentralTrack* trk;
  emcClusterContainer* emccont;
  emcClusterContent* emc;
  ReactionPlaneObject* rp;
  ReactionPlaneSngl* rpsngl;      
  RunHeader* header;
	
  THmulf* eoverp; // E/p calibration
  THmulf* dep;
  THmulf* emcdphi; // emcdphi calibration
  THmulf* emcdz; // emcdz calibration
  
  THmulf* bbczdc; // to check centrality, bbc mean & sigma
  THmulf* bbcrp; // to check reaction plane 
  TH1F* h_ntrk; // trk multiplicity
  TH1F* h_ncharged; // electron/positron multiplicity
  TH1F* h_nclust; // cluster mulitiplicity
  TH1F* h_nemc; // emc photon multiplicity
  
  // dead maps
  THmulf* map_dc;
  THmulf* map_emc;
  // THmulf* map_rich;
       
  int runnumber;
  int runselection;
  int systemselection;
        
  string OutFileName;  
};

#endif /* __MPHOTONEVENTQA_H__ */
