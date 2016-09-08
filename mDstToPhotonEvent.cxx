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
#include "dpReco.h"

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

#include "PhotonEvent.h"
#include "PHPhotonEvent.h"
//#include "Run16dAuCut.h"
#include "mDstToPhotonEvent.h"

//=============================================================
mDstToPhotonEvent::mDstToPhotonEvent(char *outfile, char* lookup_file) : 
  SubsysReco("mDstToPhotonEvent"),
  fFile(NULL),
  fTree(NULL),
  fEvent(NULL)
  // docut(NULL)
{
  fReconstruction = new dpReco(lookup_file);
  fOutFileName = outfile;
  if(fOutFileName.Length()>0)
    fTree = new TTree("T","one tree to analyze them all");
  fLookupFileName = lookup_file;
  // std::cout << "mDstToPhotonEvent::CTOR()" << std::endl;
  return;
}
//=============================================================
mDstToPhotonEvent::~mDstToPhotonEvent()
{
  if(fTree) delete fTree;
  if(fFile) delete fFile;
  // if(docut) delete docut;
}
//=============================================================
int mDstToPhotonEvent::Init(PHCompositeNode *)
{
  std::cout << "mDstToPhotonEvent::INIT()" << std::endl;
  return 0;
}
//=============================================================
int mDstToPhotonEvent::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(nodeIter.findFirst("PHCompositeNode","DST"));
  if(dstNode == NULL) {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHPhotonEvent *phevent = new PHPhotonEvent();
  PHIODataNode<PHPhotonEvent> *hitNode = new PHIODataNode<PHPhotonEvent>(phevent,"PhotonEvent","PHPhotonEvent");
  dstNode->addNode(hitNode);
  fEvent = phevent->GetEvent();
  if(fOutFileName.Length()>0) {
    fTree->Branch("PhotonEvent",fEvent);
    fTree->Reset();
  }
  // docut=new Run16dAuCut();
  //std::cout << "mDstToPhotonEvent::INITRUN()" << std::endl;
  return 0;
}
//=============================================================
int mDstToPhotonEvent::process_event(PHCompositeNode *topNode)
{
  fEvent->Clear();
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
  fEvent->SetRunNumber( header->get_RunNumber() );
  fEvent->SetBBCcharge( phg->getBbcChargeN() + phg->getBbcChargeS() );
  fEvent->SetVtxZ( phg->getBbcZVertex() );
  fEvent->SetCentrality( phg->getCentrality() );

  ReactionPlaneSngl *rpsngl;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 0, 1));
  float psi_BBCS = (rpsngl) ? rpsngl->GetPsi() : -9999;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 1, 1));
  float psi_BBCN = (rpsngl) ? rpsngl->GetPsi() : -9999;
  rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1));
  float psi_BBC = (rpsngl) ? rpsngl->GetPsi() : -9999;
  fEvent->SetPsiBBC(psi_BBC);
  fEvent->SetPsiBBCN(psi_BBCN);
  fEvent->SetPsiBBCS(psi_BBCS);

  // PHPoint vtxr = d_vtxout->get_Vertex("SVX_PRECISE_RECAL");
  // event.SetSvxZ( vtxr.getZ() );

  // EVENT CUTS
  if(cnttrk->get_npart()<2)return DISCARDEVENT;
  //trigger selection
  unsigned int trigger_scaled = _Trig_ptr->get_lvl1_trigscaled();
  unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
  unsigned int trigger_FVTXNSBBCS        = 0x00400000;
  unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
  unsigned int trigger_BBCLL1narrow      = 0x00000010;
  unsigned int accepted_triggers = 0;
  int runnumber=header->get_RunNumber();
  // --- Run16dAu200                                                                                   
  if ( runnumber >= 454774 && runnumber <= 455639 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
  // --- Run16dAu62                                                                                                        
  if ( runnumber >= 455792 && runnumber <= 456283 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
  // --- Run16dAu20                                                                                                        
  if ( runnumber >= 456652 && runnumber <= 457298 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
  // --- Run16dAu39                                                                                                        
  if ( runnumber >= 457634 && runnumber <= 458167 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
  unsigned int trigscaled_on = trigger_scaled & accepted_triggers;
  if(!trigscaled_on)return DISCARDEVENT;

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
    float Z_GLOBAL=75,EMCDZ=20,EMCDPHI=0.05,N0=1,DISP=5,CHI2_NEP0=10,PROB=0.01;//DEP0=-2,DEP1=2;
    if ( ((cnttrk->get_quality(t)==31) || (cnttrk->get_quality(t)==51) || (cnttrk->get_quality(t)==63 ))
	 && (fabs(cnttrk->get_zed(t))< Z_GLOBAL)
	 && (fabs(cnttrk->get_emcdz(t))< EMCDZ)
	 && (fabs(cnttrk->get_emcdphi(t))< EMCDPHI)
	 && (cnttrk->get_n0(t)> N0)
	 && (cnttrk->get_disp(t)<DISP)
	 && ((cnttrk->get_chi2(t)/cnttrk->get_npe0(t))<CHI2_NEP0)
	 && (cnttrk->get_prob(t)>PROB)){
	 // &&(cnttrk->get_dep(t)>DEP0)
	 // && (cnttrk->get_dep(t)< DEP1)){
    //after cut do this
      if(track.GetCharge()==1)fEvent->AddNTrack( track );
      if(track.GetCharge()==-1)fEvent->AddPTrack( track );
    }
  }
  
  // PAIRS
  CNTE *ele, *pos;
  bool atLeastOneFound = false;
  float PHIV=0.1;
  float DZ_DC=4;
  float R_CONV[2];
  R_CONV[0]=1;
  R_CONV[1]=29;
  float DPHI_CONV=0.01;//original 0.005 for AuAu
  for(int in=0; in!=fEvent->GetNEtracks(); ++in) // electrons
    for(int ip=0; ip!=fEvent->GetNPtracks(); ++ip) { // positrons
      ele = fEvent->GetNtrack(in); // get electron
      pos = fEvent->GetPtrack(ip); // get positron
      CNTDE *de = new CNTDE(); // dielectron
      fReconstruction->findIntersection(ele,pos,de,fEvent->GetVtxZ()); // FIXME vertexZ
      // PAIR CUTS
      float phiv = getPhiv(ele->GetPx(),ele->GetPy(),ele->GetPz(), pos->GetPx(), pos->GetPy(), pos->GetPz() );
      if ( fabs(phiv-TMath::Pi())>=PHIV ) continue; // phiv cut                          
      if ( fabs(ele->GetZed()-pos->GetZed())>=DZ_DC ) continue;//dz cut
      //conversion pair check
      float radius_r = de->GetR();
      float dphi_r   = de->GetEphi()-de->GetPphi();
      if ( dphi_r>TMath::Pi() )  dphi_r = 2*TMath::Pi()-dphi_r; // 1.3pi->0.7pi
      if ( dphi_r<-TMath::Pi() ) dphi_r = -2*TMath::Pi()-dphi_r; // -1.3pi->-0.7pi
      // float dtheta_r = pair.GetThetaElectron()-pair.GetThetaPositron();
      // std::cout << " VEROPAIR " << radius_r << " " << dphi_r << std::endl;
      if ( radius_r<=R_CONV[0] || radius_r>=R_CONV[1] || fabs(dphi_r)>=DPHI_CONV) continue;
      atLeastOneFound = true;
      delete de;
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
    if ( emc->chi2()< 3 ){
      fEvent->AddCluster( cluster );
    }
  }

  if(fOutFileName.Length()>0) fTree->Fill();

  return 0;
}
//=============================================================
float mDstToPhotonEvent::getPhiv(float px_e, float py_e, float pz_e, float px_p, float py_p, float pz_p)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;
  double Me=0.000510998918;
  double Me2=Me*Me;
  p1.SetX(px_e);
  p1.SetY(py_e);
  p1.SetZ(pz_e);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));
  p2.SetX(px_p);
  p2.SetY(py_p);
  p2.SetZ(pz_p);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));
  pair = p1 + p2; // pair corresponds to photon if the pair matches                                                     

  TVector3 P1, P2, Photon, z, v, u, w, wc;
  z.SetX(0);
  z.SetY(0);
  z.SetZ(1); // unit vector along z                                                                                     

  P1.SetX(px_e);
  P1.SetY(py_e);
  P1.SetZ(pz_e);

  P2.SetX(px_p);
  P2.SetY(py_p);
  P2.SetZ(pz_p);

  Photon.SetX(pair.Px());
  Photon.SetY(pair.Py());
  Photon.SetZ(pair.Pz());

  v = (P1.Cross(P2)).Unit(); // unit vector v corresponds to the unit vector of the cross product of ep pair            
  u = Photon.Unit();
  w = (u.Cross(v)).Unit();
  wc = (u.Cross(z)).Unit();

  return acos(-w.Dot(wc));
}

//=============================================================
int mDstToPhotonEvent::End(PHCompositeNode *topNode)
{
  //std::cout<<"mDstToPhotonEvent::End() ==> Terminating... "<<std::endl;
  if(fOutFileName.Length()>0) {
    std::cout << "FileName: " << fOutFileName.Data() << std::endl;
    fFile = new TFile(fOutFileName.Data(),"RECREATE");
    fTree->Write();
    fFile->Write();
    fFile->Close();
  }
  //std::cout<<"[DONE]"<<std::endl;
  return 0;
}
