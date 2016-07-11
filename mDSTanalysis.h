#ifndef __MDSTANALYSIS_H__
#define __MDSTANALYSIS_H__

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "SubsysReco.h"
#include "PhotonEvent.h"
#include "Reconstruction.h"
#include "MyCut.h"

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

class mDSTanalysis: public SubsysReco {
 public:
  mDSTanalysis(char *outfile = "output.root", char *lut = "lookup.root");
  virtual ~mDSTanalysis();
  int   Init(PHCompositeNode*);
  int   InitRun(PHCompositeNode*);
  int   process_event(PHCompositeNode*);
  int   ResetEvent(PHCompositeNode*) {return 0;};
  int   Reset(PHCompositeNode*) {return 0;};
  int   End(PHCompositeNode*);

  int   GetNodes(PHCompositeNode*);
  void  fill_trks_to_myevent();
  int   check_myevent_for_myrecopairs();
  void  fill_emcclusters_to_myevent();
  
 protected:   
  TFile *fFile;
  TTree fT;
  CentralityReco *fCentrality;
  Reconstruction fReconstruction;
  PhotonEvent fEvent;
  MyCut *fCut;
  
  PHGlobal* phg;
  TrigLvl1* _Trig_ptr;    
  PHCentralTrack* trk;
  emcClusterContainer* emccont;
  emcClusterContent* emc;
  ReactionPlaneObject* rp;
  ReactionPlaneSngl* rpsngl;
  RunHeader* header;
  
  int fRunNumber;
  int fPairCutSelection; // if ( paircutselection == 1 ) apply all pair cuts
  
  TString fOutFileName;
  TString fLookupFileName; 
}

#endif /* __RDANALYZER_H__ */



