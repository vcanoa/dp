#include <sys/resource.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include "TTree.h"
#include <TGraph.h>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Reconstruction.h"

#include "phool.h"

//Global Includes
#include "getClass.h"
#include "PHGlobal.h"
#include "recoConsts.h"
#include "Fun4AllServer.h"
#include "Fun4AllReturnCodes.h"
#include "ReactionPlaneObject.h"
#include "ReactionPlaneSngl.h"
#include "RpConst.h"
// #include "VtxOut.h"

//Trigger/Header Includes
#include "utiCentrality.h"
#include "TrigLvl1.h"
#include "TriggerHelper.h"
#include "RunHeader.h"

//Track Includes
#include "PHSnglCentralTrack.h"
#include "PHCentralTrack.h"
#include "PHCentralTrackv24.h"
#include "PHSnglCentralTrackv24.h"

//EMCal includes
#include "emcClusterContainer.h"
#include "emcClusterContent.h"


#include "mDSTanalysis.h"

//=============================================================
mDSTanalysis::mDSTanalysis(char *outfile, char* lookup_file) : 
  SubsysReco("mDSTanalysis"),
  fFile(NULL),
  fTree("T","one tree to analyze them all")
{
  fReconstruction = new Reconstruction();
  fEvent = PhotonEvent();
  fOutFileName = outfile;
  fLookupFileName = lookup_file;
  return;
}
//=============================================================
mDSTanalysis::~mDSTanalysis()
{
  if(fFile) delete fFile;
}
//=============================================================
int mDSTanalysis::Init(PHCompositeNode *)
{
  fFile = new TFile(fOutFileName.Data(),"RECREATE");
  fTree.Branch("PhotonEvent",&fEvent);
  return 0;
}
//=============================================================
int mDSTanalysis::InitRun(PHCompositeNode *topNode)
{
  fTree.Reset();
  return 0;
}
//=============================================================
int mDSTanalysis::process_event(PHCompositeNode *topNode)
{
  static int ncalls = 0;
  ncalls++;
  if ( (ncalls)%1000==0 ) std::cout << "Event process = " << ncalls << std::endl;

  PHGlobal *phg = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
  if( !phg ) return DISCARDEVENT;
  TrigLvl1 *_Trig_ptr = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if( !_Trig_ptr ) return DISCARDEVENT;
  PHCentralTrack *cnttrk = findNode::getClass<PHCentralTrack>(topNode,"PHCentralTrack");
  if( !cnttrk ) return DISCARDEVENT;
  emcClusterContainer *emcclu = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
  if( !emcclu ) return DISCARDEVENT;
  ReactionPlaneObject *rp = findNode::getClass<ReactionPlaneObject>(topNode, "ReactionPlaneObject");
  if( !rp ) return DISCARDEVENT;
  RunHeader *header = findNode::getClass<RunHeader>(topNode,"RunHeader");
  if( !header ) return DISCARDEVENT;
  
  // EVE
  fEvent.SetRunNumber( header->get_RunNumber() );
  fEvent.SetBBCcharge( phg->getBbcChargeN() + phg->getBbcChargeS() );
  fEvent.SetVtxZ( phg->getBbcZVertex() );
  fEvent.SetCentrality( phg->getCentrality() );

  ReactionPlaneSngl *rpsngl;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 0, 1));
  float psi_BBCS = (rpsngl) ? rpsngl->GetPsi() : -9999;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 1, 1));
  float psi_BBCN = (rpsngl) ? rpsngl->GetPsi() : -9999;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1));
  float psi_BBC = (rpsngl) ? rpsngl->GetPsi() : -9999;
  fEvent.SetPsiBBC(psi_BBC);
  fEvent.SetPsiBBCN(psi_BBCN);
  fEvent.SetPsiBBCS(psi_BBCS);

  // PHPoint vtxr = d_vtxout->get_Vertex("SVX_PRECISE_RECAL");
  // event.SetSvxZ( vtxr.getZ() );

  // EVENT CUTS
  
  // TRACKS
  for (int t=0; t!=int( cnttrk->get_npart() ); ++t) {
    CNTE track;
    track.SetDCarm( cnttrk->get_dcarm(t) );
    track.SetDCside( cnttrk->get_dcside(t) );
    track.SetPx( cnttrk->get_px(t) );
    track.SetPy( cnttrk->get_py(t) );
    track.SetPz( cnttrk->get_pz(t) );
    track.SetPhi( cnttrk->get_phi(t) );
    track.SetPhi0( cnttrk->get_phi0(t) );
    track.SetTheta0( cnttrk->get_the0(t) );
    track.SetZed( cnttrk->get_zed(t) );
    track.SetAlpha( cnttrk->get_alpha(t) );
    track.SetCharge( cnttrk->get_charge(t) );
    track.SetEMCid( cnttrk->get_emcid(t) );
    track.SetEcore( cnttrk->get_ecore(t) );
    track.SetDisp( cnttrk->get_disp(t) );
    track.SetN0( cnttrk->get_n0(t) );
    track.SetChi2( cnttrk->get_chi2(t) );
    track.SetNpe0( cnttrk->get_npe0(t) );
    track.SetDep( cnttrk->get_dep(t) );
    track.SetProb( cnttrk->get_prob(t) );
    track.SetEMCdz( cnttrk->get_emcdz(t) );
    track.SetEMCdphi( cnttrk->get_emcdphi(t) );
    track.SetEMCx( cnttrk->get_pemcx(t) );
    track.SetEMCy( cnttrk->get_pemcy(t) );
    track.SetEMCz( cnttrk->get_pemcz(t) );
    track.SetCenterPhi( cnttrk->get_center_phi(t) );
    track.SetCenterZ( cnttrk->get_center_z(t) );
    track.SetPPC1x( cnttrk->get_ppc1y(t) );
    track.SetPPC1y( cnttrk->get_ppc1x(t) );
    track.SetPPC1z( cnttrk->get_ppc1z(t) );
    track.SetEsect( cnttrk->get_sect(t) );
    track.SetYsect( cnttrk->get_ysect(t) );
    track.SetZsect( cnttrk->get_zsect(t) );

    //TRACK CUTS

    if(track.GetCharge()==1)
      fEvent.AddNTrack( track );
    if(track.GetCharge()==-1)
      fEvent.AddPTrack( track );
  }

  // PAIRS
  CNTE *ele, *pos;
  bool atLeastOneFound = false;
  for(int in=0; in!=fEvent.GetNEtracks(); ++in) // electrons
    for(int ip=0; ip!=fEvent.GetNPtracks(); ++ip) { // positrons
      ele = fEvent.GetNtrack(in); // get electron
      pos = fEvent.GetPtrack(ip); // get positron
      CNTDE de; // dielectron
      //fReconstruction->findIntersection(ele,pos,de,fEvent.GetVtxZ()); // FIXME vertexZ

      // PAIR CUTS

      atLeastOneFound = true;
    }
  if(!atLeastOneFound) return DISCARDEVENT;

  // CLUSTERS
  for(int c=0; c!=int( emcclu->size() ); ++c) {
    emcClusterContent *emc = emcclu->getCluster(c);
    EMCC cluster;
    cluster.SetArm( emc->arm() );
    cluster.SetID( emc->id() );
    cluster.SetX( emc->x() );
    cluster.SetY( emc->y() );
    cluster.SetZ( emc->z() );
    cluster.SetEcore( emc->ecore() );
    cluster.SetProb( emc->prob_photon() );
    cluster.SetEmcdz( emc->emctrkdz() );
    cluster.SetEmcdphi( emc->emctrkdphi() );

    // CLUSTER CUTS

    fEvent.AddCluster( cluster );
  }

  fTree.Fill();
  fEvent.Clear();

  return 0;
}
//=============================================================
int mDSTanalysis::End(PHCompositeNode *topNode)
{
  std::cout<<"mDSTanlysis::End() ==> Terminating... "<<std::endl;
  fTree.Write();
  fFile->Write();
  fFile->Close();
  std::cout<<"[DONE]"<<std::endl;
  return 0;
}
