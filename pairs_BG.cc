#include <math.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMatrixD.h"
#include "TNtuple.h"
#include "PhotonEvent.h"
#include "dpReco.h"
#include "CNTE.h"
#include "CNTDE.h"

using namespace std;

static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;
const int cutbin = 5;
TH2F *hPhMassPt[cutbin];
TH3F *hPhMassPtPhiv;
TH3F *hPhMassPtDzed;
TH3F *hPhMassPtRad;
TH3F *hPhMassPtDthe;
TH2F *hPhoton_PT1_PT2;
TH2F *hPhoton_PT1_PT3;
TH2F *hPhotonMass1Pt1;
TH2F *hPhotonMass2Pt1;
TH2F *hPhotonMass3Pt1;
TH1F *hTrkEMCdz[2];
TH1F *hTrkEMCdphi[2];
TH1F *hTrkZed[2];
TH1F *hTrkN0[2];
TList *myList;

//=======================================
void Histograms() {
  for(int i=0; i<cutbin; ++i) hPhMassPt[i] = new TH2F(Form("hPhMassPt_cut=%d",i),";m_{ee} [GeV];p_{T}^{ee}",100,0,1, 100,0,10);
  hPhMassPtPhiv = new TH3F("hPhMassPtPhiv",";MASS_VTX [GeV];PT_VTX [rad];PhiV [rad]",  100,0,1, 100,0,10, 100,0,TMath::Pi());
  hPhMassPtDzed = new TH3F("hPhMassPtDzed",";MASS_VTX [GeV];PT_VTX [rad];dZed [cm]",   100,0,1, 100,0,10, 100,0,10);
  hPhMassPtRad  = new TH3F("hPhMassPtRad", ";MASS_VTX [GeV];PT_VTX [rad];Rad [cm]",    100,0,1, 100,0,10, 100,0,30);
  hPhMassPtDthe = new TH3F("hPhMassPtDthe",";MASS_VTX [GeV];PT_VTX [rad];dTheta [Rad]",100,0,1, 100,0,10, 100,0,0.05);
  hPhoton_PT1_PT2 = new TH2F("hPhoton_PT1_PT2",";PT_VTX [GeV];PT_PROJ [GeV]",  100,0,10, 100,0,10);
  hPhoton_PT1_PT3 = new TH2F("hPhoton_PT1_PT2",";PT_VTX [GeV];PT_CONV [GeV]",  100,0,10, 100,0,10);
  hPhotonMass1Pt1 = new TH2F("hPhotonMass1Pt1",";MASS_VTX [GeV];PT_VTX [GeV]", 100,0,1,  100,0,10);
  hPhotonMass2Pt1 = new TH2F("hPhotonMass2Pt1",";MASS_PROJ [GeV];PT_VTX [GeV]",100,0,1,  100,0,10);
  hPhotonMass3Pt1 = new TH2F("hPhotonMass3Pt1",";MASS_CONV [GeV];PT_VTX [GeV]",100,0,1,  100,0,10);
  hTrkEMCdz[0] = new TH1F("hTrkEMCdz0","Negative;DZED [cm]",100,0,10);
  hTrkEMCdz[1] = new TH1F("hTrkEMCdz1","Positive;DZED [cm]",100,0,10);
  hTrkEMCdphi[0] = new TH1F("hTrkEMCdphi0","Negative;DPHI [rad]",100,0,0.1);
  hTrkEMCdphi[1] = new TH1F("hTrkEMCdphi1","Positive;DPHI [rad]",100,0,0.1);
  hTrkZed[0] = new TH1F("hTrkZed0","Negative;ZED [cm]",100,0,100);
  hTrkZed[1] = new TH1F("hTrkZed1","Positive;ZED [cm]",100,0,100);
  hTrkN0[0] = new TH1F("hTrkNZERO0","Negative;N0",100,0,100);
  hTrkN0[1] = new TH1F("hTrkNZERO1","Positive;N0",100,0,100);
  for(int i=0; i!=2; ++i) {
    myList->Add( hTrkEMCdz[i] );
    myList->Add( hTrkEMCdphi[i] );
    myList->Add( hTrkZed[i] );
    myList->Add( hTrkN0[i] );
  }
  for(int i=0; i<cutbin; ++i) myList->Add( hPhMassPt[i] );
  myList->Add( hPhMassPtPhiv );
  myList->Add( hPhMassPtDzed );
  myList->Add( hPhMassPtRad );
  myList->Add( hPhMassPtDthe );
  myList->Add( hPhoton_PT1_PT2 );
  myList->Add( hPhoton_PT1_PT3 );
  myList->Add( hPhotonMass1Pt1 );
  myList->Add( hPhotonMass2Pt1 );
  myList->Add( hPhotonMass3Pt1 );
}
//=======================================
void getPtAndMass(float& pt, float& ms, float px1, float py1, float pz1, float px2, float py2, float pz2) {
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;
  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));
  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));
  pair = p1 + p2;
  pt = pair.Pt();
  ms = pair.M();
}
//=======================================
float getPhiv(CNTE *neg, CNTE *pos) {
  float px_e = neg->GetPx();
  float py_e = neg->GetPy();
  float pz_e = neg->GetPz();
  float px_p = pos->GetPx();
  float py_p = pos->GetPy();
  float pz_p = pos->GetPz();
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;
  p1.SetX(px_e);
  p1.SetY(py_e);
  p1.SetZ(pz_e);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));
  p2.SetX(px_p);
  p2.SetY(py_p);
  p2.SetZ(pz_p);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));
  pair = p1 + p2;
  TVector3 P1, P2, Photon, z, v, u, w, wc;
  z.SetX(0);
  z.SetY(0);
  z.SetZ(1);
  P1.SetX(px_e);
  P1.SetY(py_e);
  P1.SetZ(pz_e);
  P2.SetX(px_p);
  P2.SetY(py_p);
  P2.SetZ(pz_p);
  Photon.SetX(pair.Px());
  Photon.SetY(pair.Py());
  Photon.SetZ(pair.Pz());
  v = (P1.Cross(P2)).Unit();
  u = Photon.Unit();
  w = (u.Cross(v)).Unit();
  wc = (u.Cross(z)).Unit();
  float phiv = acos(-w.Dot(wc));
  return phiv;
}
//=======================================
bool passtrkcuts(CNTE *trk) {
  bool pass = true;
  if( TMath::Abs(trk->GetEMCdz()) > 10 ) pass=false;
  if( TMath::Abs(trk->GetEMCdphi()) > 0.025 ) pass=false;
  if( TMath::Abs(trk->GetZed()) > 75 ) pass=false;
  if(trk->GetN0()<=3 ) pass=false;
  return pass;
}
//=======================================
void qatrack(CNTE *trk) {
  int lvl = trk->GetCharge()<0?0:1;
  hTrkEMCdz[lvl]->Fill( trk->GetEMCdz() );
  hTrkEMCdphi[lvl]->Fill( trk->GetEMCdphi() );
  hTrkZed[lvl]->Fill( trk->GetZed() );
  hTrkN0[lvl]->Fill( trk->GetN0() );
}
//=======================================
void readgoodtracks(PhotonEvent *event, int& npos, int& nneg, int pos[20], int neg[20], bool qa=false) {
  npos = 0;
  nneg = 0;
  CNTE *track1;
  for(int k=0; k<event->GetNPtracks(); ++k) {
    track1 = event->GetPtrack(k);
    if(qa) qatrack( track1 );
    if(!passtrkcuts(track1)) continue;
    pos[npos++] = k;
    if(npos==20) {
      printf("more than 20 tracks in event... ignoring the rest!!!\n");
      break;
    }
  }
  for(int k=0; k<event->GetNEtracks(); ++k) {
    track1 = event->GetNtrack(k);
    if(qa) qatrack( track1 );
    if(!passtrkcuts(track1)) continue;
    neg[nneg++] = k;
    if(nneg==20) {
      printf("more than 20 tracks in event... ignoring the rest!!!\n");
      break;
    }
  }
}
//=======================================
void BuildPair(CNTE *track1, CNTE *track2, dpReco *reco, float zVtx) {
  //============
  // PAIRS cuts
  float DTHETA = 0.1;
  float RCONV[2] = {9,23};
  //============
  // QUICK PAIR CUTS
  float ZEDCUT = 4;
  float PHIVCUT = 0.1;
  //===========================

  // BUILDING PAIR
  CNTDE mypair;
  CNTE mytrk1(track1); //make obj from pointer to obj
  CNTE mytrk2(track2); //make obj from pointer to obj
  float px_1 = mytrk1.GetPx();
  float py_1 = mytrk1.GetPx();
  float pz_1 = mytrk1.GetPx();
  float px_2 = mytrk1.GetPx();
  float py_2 = mytrk1.GetPx();
  float pz_2 = mytrk1.GetPx();
  float PT1_1 = TMath::Sqrt( px_1*px_1 + py_1*py_1 );
  float PT1_2 = TMath::Sqrt( px_2*px_2 + py_2*py_2 );
  reco->findIntersection(&mytrk1, &mytrk2, &mypair, zVtx);
  float rad = mypair.GetR();
  float phi = (mypair.GetEphi()+mypair.GetPphi())/2;
  float the = (mypair.GetEtheta()+mypair.GetPtheta())/2;
  float dthe = mypair.GetEtheta()-mypair.GetPtheta();
  dthe = TMath::Abs(dthe);
  TVector3 mom1 = reco->findMomentum(&mytrk1, rad, phi, the, zVtx);
  TVector3 mom2 = reco->findMomentum(&mytrk2, rad, phi, the, zVtx);
  float PT2_1 = mom1.Pt();
  float PT2_2 = mom2.Pt();
  bool pass_phiv = false;
  bool pass_dzed = false;
  bool pass_rad  = false;
  bool pass_dthe = false;
  // PHIV
  float phiv = getPhiv(track2,track1); // neg pos (convention)
  phiv = TMath::Abs(phiv-TMath::Pi());
  if( phiv<=PHIVCUT ) pass_phiv = true;
  // dZED
  float dzed = track2->GetZed() - track1->GetZed(); // neg - pos (convention)
  dzed = TMath::Abs(dzed);
  if( dzed<=ZEDCUT ) pass_dzed = true;
  // RAD
  if( (rad>=RCONV[0]) && (rad<=RCONV[1]) ) pass_rad = true;
  // DTHETA
  if( dthe<=DTHETA ) pass_dthe = true;
  
  // COMPUTING PAIR OBSERVABLES
  float phPt1, phPt2, phPt3;
  float phMass1, phMass2, phMass3;
  float cp = TMath::Cos(phi);
  float sp = TMath::Sin(phi);
  getPtAndMass(phPt1,phMass1, px_2,     py_2,     pz_2, px_1,     py_1,     pz_1);
  getPtAndMass(phPt2,phMass2, PT1_2*cp, PT1_2*sp, pz_2, PT1_1*cp, PT1_1*sp, pz_1);
  getPtAndMass(phPt3,phMass3, PT2_2*cp, PT2_2*sp, pz_2, PT2_1*cp, PT2_1*sp, pz_1);
  // FILLING CONTROL HISTOGRAMS
  hPhMassPtPhiv->Fill(phMass1, phPt1, phiv);
  hPhMassPtDzed->Fill(phMass1, phPt1, dzed);
  hPhMassPtRad->Fill(phMass1, phPt1, rad);
  hPhMassPtDthe->Fill(phMass1, phPt1, dthe);
  if( pass_phiv ) hPhMassPt[0]->Fill(phMass1, phPt1);
  if( pass_dzed ) hPhMassPt[1]->Fill(phMass1, phPt1);
  if( pass_rad ) hPhMassPt[2]->Fill(phMass1, phPt1);
  if( pass_dthe ) hPhMassPt[3]->Fill(phMass1, phPt1);
  // FILLING RESULT HISTOGRAMS
  if( !pass_phiv || !pass_dzed || !pass_rad || !pass_dthe ) return;
  hPhoton_PT1_PT2->Fill(phPt1,phPt2);
  hPhoton_PT1_PT3->Fill(phPt1,phPt3);
  hPhotonMass1Pt1->Fill(phMass1, phPt1);
  hPhotonMass2Pt1->Fill(phMass2, phPt1);
  hPhotonMass3Pt1->Fill(phMass3, phPt1);
}
//=======================================
int main(int nvar, char** vars) {
  //============
  // CONNECTION
  int run=atoi(vars[1]); // parameter 1
  int seg=atoi(vars[2]); // parameter 2
  const char* inFile =Form("outputPhoton/dAu200GeV_%010d-%04d.root",run,seg);  
  dpReco reco( "lookup_3D.root" );
  TFile* f1 = new TFile(inFile,"READ");
  TTree* T = (TTree*)f1->Get("T");
  TBranch* br = T->GetBranch("PhotonEvent");
  PhotonEvent *event=0;
  PhotonEvent *preevent=0;
  br->SetAddress(&event);
  //============
  
  int nevt = T->GetEntries();
  for(int i1=0; i1<nevt; ++i1) {
    if(i1%1000 == 0) cout << "Event:  " << i1 << "/" << nevt << endl;
    br->GetEntry(i1);
    float zVtx = event->GetVtxZ(); 
    int ntrk = event->GetNPtracks() + event->GetNEtracks();
    if(ntrk==0) continue; // discard event, but dont touch previous event

    // ==> Filtering good tracks in this event <==
    int npos=0, pos[20];
    int nneg=0, neg[20];
    readgoodtracks(event,npos,nneg,pos,neg, true);
    if((npos+nneg)==0) continue; // discard event, but dont touch previous event
    if(!preevent) { // setting up previous event for the first time
      preevent = event;
      continue;
    }

    // ==> Filtering good tracks in previous event <==
    int nprepos=0, prepos[20];
    int npreneg=0, preneg[20];
    readgoodtracks(preevent,nprepos,npreneg,prepos,preneg);
    
    // ==> Making pairs <==
    CNTE *track1, *track2;
    for(int ip=0; ip<npos; ++ip) {
      track1 = event->GetPtrack( pos[ip] ); // getting one good positive track from this event
      for(int jn=0; jn<npreneg; ++jn) {
	track2 = preevent->GetNtrack( preneg[jn] ); // getting one good negative track from previous event
	if( track1->GetDCarm()!=track2->GetDCarm() ) continue;
	BuildPair(track1,track2,&reco,zVtx);
      }
    }
    for(int in=0; in<nneg; ++in) {
      track2 = event->GetNtrack( neg[in] ); // getting one good negative track from this event
      for(int jp=0; jp<nprepos; ++jp) {
	track1 = preevent->GetPtrack( prepos[jp] ); // getting one good positive track from previous event
	if( track1->GetDCarm()!=track2->GetDCarm() ) continue;
	BuildPair(track1,track2,&reco,zVtx);
      }
    }
  }
  f1->Close();

  TFile* f2 = new TFile(Form("output_FG/dAu200_mixingbg_%d_%d.root",run,seg),"RECREATE");
  for(int i=0; i!=myList->GetEntries(); ++i) myList->At(i)->Write();
  f2->Close();

  return 0;
}
