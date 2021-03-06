#ifndef __MDSTTOPHOTONEVENT_H__
#define __MDSTTOPHOTONEVENT_H__

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "SubsysReco.h"

#include "PhotonEvent.h"
#include "dpReco.h"

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

class mDstToPhotonEvent : public SubsysReco {
 public:
  mDstToPhotonEvent(char *outfile, char *lut);
  virtual ~mDstToPhotonEvent();
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
  TTree *fTree;
 
  PhotonEvent *fEvent;
  dpReco *fReconstruction;
  
  TString fOutFileName;
  TString fLookupFileName; 
};

#endif /* __MDSTTOPHOTONEVENT_H__ */
