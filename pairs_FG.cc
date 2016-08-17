#include <math.h>
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TNtuple.h"
#include "PhotonEvent.h"
#include "dpReco.h"
#include "CNTE.h"
#include "CNTDE.h"

using namespace std;

static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;
static const double c = 299792458;

float getPT(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
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

  pair = p1 + p2;//pair corresponds to photon if the pair matches

  return pair.Pt();
}

float getMass(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
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
 
  return pair.M();
}

float getPi0Mass(float px1, float py1, float pz1, float px2, float py2, float pz2, float px, float py, float pz)
{// px, py, pz correspond to emcal photon cluster 
  TLorentzVector p1, p2;
  TLorentzVector convPhoton, emcPhoton;
  TLorentzVector pi0;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  convPhoton = p1 + p2;

  emcPhoton.SetX(px);
  emcPhoton.SetY(py);
  emcPhoton.SetZ(pz);
  emcPhoton.SetE(emcPhoton.P());

  pi0 = convPhoton+emcPhoton;

  return pi0.M();
}

float getPi0Pt(float px1, float py1, float pz1, float px2, float py2, float pz2, float px, float py, float pz)
{// px, py, pz correspond to emcal photon cluster 
  TLorentzVector p1, p2;
  TLorentzVector convPhoton, emcPhoton;
  TLorentzVector pi0;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  convPhoton = p1 + p2;

  emcPhoton.SetX(px);
  emcPhoton.SetY(py);
  emcPhoton.SetZ(pz);
  emcPhoton.SetE(emcPhoton.P());

  pi0 = convPhoton+emcPhoton;

  return pi0.Pt();
}

float getConvertedPhotonE(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
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

  return pair.E();
}

float getDcenter(float phi1, float z1, float phi2, float z2)
{
  float dcenter_phi_sigma = 0.01, dcenter_phi_offset = 0.;
  float dcenter_z_sigma = 3.6, dcenter_z_offset = 0.;

  float dcenter_z = (z1-z2-dcenter_z_offset)/dcenter_z_sigma;
  float dcenter_phi = (phi1-phi2-dcenter_phi_offset)/dcenter_phi_sigma;

  return sqrt(dcenter_phi*dcenter_phi+dcenter_z*dcenter_z);
}

float getPhiv(float px_e, float py_e, float pz_e, float px_p, float py_p, float pz_p)
{
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

  float phiv = acos(-w.Dot(wc));
  return phiv;
}

int main()
{
  const char* inFile = "dAu200GeV_0000454774-9000.root";
  int job = 0;

  //gSystem->Load("libDPRun16dAu.so");
  dpReco reco( "lookup_3D.root" );

  TFile* f1 = new TFile(inFile,"READ");
  if(!(f1)) {
    cout<<"can't find input file..."<<endl;
    exit(1);
  }
  TTree* T = (TTree*)f1->Get("T");
  TBranch* br = T->GetBranch("PhotonEvent");
  PhotonEvent* event=0;
  br->SetAddress(&event);
  
  const int centbin = 5;
  float centlow[centbin] = {0, 20, 40, 60, 80};
  float centup[centbin]  = {20, 40, 60, 80, 100};

  //======================================================
  //                    diagnostic
  //======================================================
  TH1F* fgPhiDC = new TH1F("fgPhiDC","#phi_{DC} distribution;#phi_{DC} (rad)",500,-0.5*TMath::Pi(),1.5*TMath::Pi());

  //======================================================
  //           fg & bg of the cutting variables
  //======================================================
  TH2F* fgPhiv[centbin] = {0};
  TH2F* fgRadius[centbin] = {0};
  TH2F* fgPhi[centbin] = {0};
  TH2F* fgTheta[centbin] = {0};
  TH2F* fgDz[centbin] = {0};
  TH2F* fgDtheta[centbin] = {0};
  // TH2F* fgDphi[centbin] = {0};
  for (int i = 0; i < centbin; ++i)
  {
    fgPhiv[i] = new TH2F(Form("fgPhiv[%d]", i),"#phi_{V} distribution; #phi_{V} (rad); p_{T}^{ee}",1000,0,TMath::Pi(),1000,0,10);
    fgRadius[i] = new TH2F(Form("fgRadius[%d]", i),"r_{conv} distribution; r_{conv} (cm); p_{T}^{ee}",1000,0,30,1000,0,10);
    fgPhi[i] = new TH2F(Form("fgPhi[%d]", i),"#phi_{conv} distribution; #phi_{conv} (rad); p_{T}^{ee}",1000,-0.5*TMath::Pi(),1.5*TMath::Pi(),1000,0,10);
    fgTheta[i] = new TH2F(Form("fgTheta[%d]", i),"#theta_{conv} distribution; #theta_{conv} (rad); p_{T}^{ee}",1000,0,TMath::Pi(),1000,0,10);
    fgDz[i] = new TH2F(Form("fgDz[%d]", i),"dz_{DC} distribution; dz_{DC} (cm); p_{T}^{ee}",100,0,20,1000,0,10);
    fgDtheta[i] = new TH2F(Form("fgDtheta[%d]", i),"#Delta#theta_{conv} distribution; #Delta#theta_{conv} (rad); p_{T}^{ee}",1000,0,0.05,1000,0,10);
    // fgDphi[i] = new TH2F(Form("fgDphi[%d]", i),"#Delta#phi_{conv} distribution; #Delta#phi_{conv} (rad); p_{T}^{ee}",1000,0,0.05,1000,0,10);
  }

	//======================================================
 	//         mass distribution for converted photon
  //======================================================
  TH1F* fgStdMass = new TH1F("fgStdMass","m_{ee} distribution;m_{ee} (GeV)",500,0,0.5);
  TH1F* fgCStdMass = new TH1F("fgCStdMass","m_{ee} distribution;m_{ee} (GeV)",500,0,0.05); // CStd means corrected mass using the reconstructed conversion angle
  TH1F* fgRecoMass = new TH1F("fgRecoMass","m_{ee} distribution;m_{ee} (GeV)",500,0,0.05);
  
  //======================================================
  //          pt distribution for converted photon
  //======================================================
  TH1F* fgPhotonPt = new TH1F("fgPhotonPt","p_{T}^{#gamma} distribution;p_{T}^{#gamma} (GeV)",500,0,10);  
  
  //======================================================
  //               integrated over all pt
  //======================================================
  const int cutbin = 8;
  TH1F* fgmassCutMBAllpt[cutbin] = {0};
  for(int i = 0; i < cutbin; i++)
    {
      fgmassCutMBAllpt[i] = new TH1F(Form("fgmassCutMBAllpt[%d]", i),"m_{ee} distribution ;m_{ee} (GeV)",500,0,1);
    }
  TH1F* fgmassCutAllpt[cutbin][centbin] = {0};
  for(int i = 0; i < cutbin; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassCutAllpt[i][j] = new TH1F(Form("fgmassCutAllpt[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV)",500,0,1);
	}
    }
  //======================================================
  //             normalizing factor check
  //======================================================
  TH2F* fgmassCut[cutbin][centbin] = {0};
  for(int i = 0; i < cutbin; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassCut[i][j] = new TH2F(Form("fgmassCut[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  //======================================================
  //             eff, BG rej, BG/FG check
  //======================================================
  const int phivbin = 14; // binw=0.025rad
  TH2F* fgmassPhiv[phivbin][centbin] = {0};
  for(int i = 0; i < phivbin; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassPhiv[i][j] = new TH2F(Form("fgmassPhiv[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int dzbin = 20; // binw=1cm
  TH2F* fgmassDz[dzbin][centbin] = {0};
  for(int i = 0; i < dzbin; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassDz[i][j] = new TH2F(Form("fgmassDz[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int rbin1 = 20; // binw=0.5cm
  TH2F* fgmassRadius1[rbin1][centbin] = {0};
  for(int i = 0; i < rbin1; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassRadius1[i][j] = new TH2F(Form("fgmassRadius1[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int rbin2 = 32;
  TH2F* fgmassRadius2[rbin2][centbin] = {0};
  for(int i = 0; i < rbin2; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassRadius2[i][j] = new TH2F(Form("fgmassRadius2[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int phibin1 = 14;
  TH2F* fgmassPhi1[phibin1][centbin] = {0};
  for(int i = 0; i < phibin1; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassPhi1[i][j] = new TH2F(Form("fgmassPhi1[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int phibin2 = 16;
  TH2F* fgmassPhi2[phibin2][centbin] = {0};
  for(int i = 0; i < phibin2; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassPhi2[i][j] = new TH2F(Form("fgmassPhi2[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int phibin3 = 8;
  TH2F* fgmassPhi3[phibin3][centbin] = {0};
  for(int i = 0; i < phibin3; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassPhi3[i][j] = new TH2F(Form("fgmassPhi3[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int phibin4 = 20;
  TH2F* fgmassPhi4[phibin4][centbin] = {0};
  for(int i = 0; i < phibin4; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassPhi4[i][j] = new TH2F(Form("fgmassPhi4[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  const int dthetabin = 25; // binw=0.001rad
  TH2F* fgmassDtheta[dthetabin][centbin] = {0};
  for(int i = 0; i < dthetabin; i++)
    {
      for (int j = 0; j < centbin; ++j)
	{
	  fgmassDtheta[i][j] = new TH2F(Form("fgmassDtheta[%d][%d]", i, j),"m_{ee} distribution ;m_{ee} (GeV); p_{T}^{ee}",500,0,1,1000,0,10);
	}
    }
  
  int nevt = T->GetEntries();
  
  float bbc = 0;
  int ntrkp = 0;
  int ntrkn = 0;
  float zVtx = 0;
  float cent = 0;
  
  int arm1=0, arm2=0;
  int trackid1=0, trackid2=0;
  int charge1, charge2;
  float alpha1, alpha2;
  float phi1, phi2, zed1, zed2;// position of DC hits
  float px1, py1, pz1, px2, py2, pz2, pT1, pT2;
  float crkz1, crkphi1, crkz2, crkphi2;
  
  float dcenter; // ghost cut
  
  float ZEDCUT = 4;
  float DTHETA = 0.01;
  float PHIVCUT = 0.1;
  float DCENTERCUT = 10;
  
  int nevt_ghost = 0;
  
  for (int i1 = 0; i1 < nevt; ++i1)
    {//event from Pool1
      cout<<"I am reading event= "<<i1<<endl;
      if(i1%1000 == 0) cout << "Event:  " << i1 << "/" << nevt << endl;
      //event->Clear();
      cout<<"after clear"<<endl;
      br->GetEntry(i1);
      cout<<"after get entry "<<i1<<endl;
      ntrkp = event->GetNPtracks();
      ntrkn = event->GetNEtracks();
      bbc = event->GetBBCcharge();
      zVtx = event->GetVtxZ(); 
      cent = event->GetCentrality(); 
      
      if ( ntrkp==0 || ntrkn==0) continue;
      // event selection for centrality
      int icent_ee = -999;
      for (int k = 0; k < centbin; ++k)
	{ cout<<"event selection for centrality ="<<k<<endl;
	  if ( (cent>=centlow[k]) && (cent<centup[k]) )  icent_ee = k;
	}
      if ( icent_ee<0 ) continue;
      //=========================================
      //        trk selection for evt
      //=========================================
      int ghost[20], realp[20];
      int allp[20];
      int alln[20];
      int nghost = 0, nrealp = 0, nrealn =0, nallp = 0, nalln=0;
      //e+
      for (int k1 = 0; k1 < ntrkp; ++k1)
	{
	  cout<<"loop positive tracks="<<k1<<endl;
	  CNTE*  track_p=event->GetPtrack(k1);
	  if ((track_p->GetEMCdz()<=-10) || (track_p->GetEMCdz()>=10)) continue;
	  if  ((track_p->GetEMCdphi()<=-0.025) || (track_p->GetEMCdphi()>=0.025)) continue;
	  if  (track_p->GetN0()<=3 ) continue;
	  if ((track_p->GetDep()<=-1.5) || (track_p->GetDep()>=2)) continue;
	  if ( fabs((track_p->GetZed())>=75 )) continue;
	  allp[nallp] = k1;
	  ++nallp;        
	  int ghost_flag = 0;
	  //e-
	  for (int k2 = 0; k2 < ntrkn; ++k2)
	    {
	      cout<<"loop negative tracks="<<k2<<endl;
	      CNTE*  track_n=event->GetNtrack(k2);
	      if ((track_n->GetEMCdz()<=-10) || (track_n->GetEMCdz()>=10)) continue;
	      if  ((track_n->GetEMCdphi()<=-0.025) || (track_n->GetEMCdphi()>=0.025)) continue;
	      if  (track_n->GetN0()<=3 ) continue;
	      if ((track_n->GetDep()<=-1.5) || (track_n->GetDep()>=2)) continue;
	      if ( fabs((track_n->GetZed())>=75 )) continue;
	      if ( getDcenter(track_p->GetCenterPhi(), track_p->GetCenterZ(), track_n->GetCenterPhi(), track_n->GetCenterZ()) < DCENTERCUT ) { ghost_flag = 1; break; }
	    }
	  if ( ghost_flag ) { ghost[nghost] = k1; ++nghost; }
	  else { realp[nrealp] = k1; ++nrealp; }
	} 
      if (nghost>0) nevt_ghost++; // # of evts that contains ghost trks
      
      //======================================================
      //           real trks X real trks/ghost trks
    //======================================================
      for (int ireal1 = 0; ireal1 < nrealp; ++ireal1)
	{
	  cout<<"loop real positive tracks= "<<ireal1<<endl;
	  CNTE*  track_p=event->GetPtrack(ireal1);
	  px1     = track_p->GetPx();
	  py1     = track_p->GetPy();
	  pz1     = track_p->GetPz();
	  pT1     = TMath::Hypot(px1, py1); // pT from std reconctruction
	  
	  phi1    = track_p->GetPhi();
	  zed1    = track_p->GetZed();
	  alpha1  = track_p->GetAlpha();
	  
	  arm1     = track_p->GetDCarm();
	  trackid1 = track_p->GetEMCid();   
	  charge1  = track_p->GetCharge();
	  
	  crkz1    = track_p->GetCenterZ();
	  crkphi1  = track_p->GetCenterPhi();
	  
	  CNTE mytrk1;
	  mytrk1.SetPhi(phi1);
	  mytrk1.SetZed(zed1);
	  mytrk1.SetAlpha(alpha1);
	  mytrk1.SetPx(px1); // standard reconstruction
	  mytrk1.SetPy(py1);
	  mytrk1.SetPz(pz1);
	  
	  fgPhiDC->Fill(phi1);
	  
	  for (int ireal2 = 0; ireal2 < nrealp; ++ireal2)
	    {
	      cout<<"loop real negative tracks= "<<ireal2<<endl;
	      CNTE*  track_n=event->GetNtrack(ireal2);
	      px2     = track_n->GetPx();
	      py2     = track_n->GetPy();
	      pz2     = track_n->GetPz();
	      pT2     = TMath::Hypot(px2, py2);
	      
	      phi2    = track_n->GetPhi();
	      zed2    = track_n->GetZed();
	      alpha2  = track_n->GetAlpha();
	      
	      arm2     = track_n->GetDCarm();
	      trackid2 = track_n->GetEMCid();
	      charge2  = track_n->GetCharge();
	      
	      crkz2    = track_n->GetCenterZ();
	      crkphi2  = track_n->GetCenterPhi();
	      
	      CNTE mytrk2;
	      mytrk2.SetPhi(phi2);
	      mytrk2.SetZed(zed2);
	      mytrk2.SetAlpha(alpha2);
	      mytrk2.SetPx(px2); // standard reconstruction
	      mytrk2.SetPy(py2);
	      mytrk2.SetPz(pz2);
	      
	      if ( arm1!=arm2 ) continue; // same arm requirement
	      if ( charge1==charge2 ) continue; // opposite charge requirement
	      
	      float phiv = 0, dzed = 0; // not sure about the initialization
	      if ( charge1==1 )//e-e+
		{//1 is electron, 2 is positron
		  phiv = getPhiv(px1, py1, pz1, px2, py2, pz2);
		  dzed = zed1-zed2;
		}
	      if ( charge1==-1 )//e+e-
		{//1 is positron, 2 is electron
		  phiv = getPhiv(px2, py2, pz2, px1, py1, pz1);
		  dzed = zed2-zed1;
		}
	      float massvtx  = getMass(px1, py1, pz1, px2, py2, pz2);
	      float phoT_std = getPT(px1, py1, pz1, px2, py2, pz2);
        
        //=====================================================
        //                      no cut 
        //=====================================================
        fgmassCutMBAllpt[0]->Fill(massvtx);
        fgmassCutAllpt[0][icent_ee]->Fill(massvtx);
        fgmassCut[0][icent_ee]->Fill(massvtx, phoT_std);

        int phiv_cut = 0, dz_cut = 0;
        if ( TMath::Abs(phiv-TMath::Pi())<=PHIVCUT ) phiv_cut = 1;
        if ( TMath::Abs(dzed)<=ZEDCUT ) dz_cut = 1;

        // if ( phiv_cut && dtheta_cut && r_cut && phi_cut ) std::cout<<"reconstructed conversion point"<<radius_r<<" "<<phi_r<<" "<<theta_r<<std::endl;

        //=====================================================
        //                     phiv cut
        //=====================================================
        fgPhiv[icent_ee]->Fill(phiv, phoT_std);
        if ( phiv_cut )
        {
          fgmassCutMBAllpt[1]->Fill(massvtx);
          fgmassCutAllpt[1][icent_ee]->Fill(massvtx);
          fgmassCut[1][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < phivbin; ++icut)
        {// phiv
          double cutphiv = 0.025*(icut+1);
          if (phiv>TMath::Pi()-cutphiv) fgmassPhiv[icut][icent_ee]->Fill(massvtx, phoT_std);
        }

        //=====================================================
        //                      dz cut
        //=====================================================
        fgDz[icent_ee]->Fill(dzed, phoT_std);
        if ( dz_cut )
        {
          fgmassCutMBAllpt[2]->Fill(massvtx);
          fgmassCutAllpt[2][icent_ee]->Fill(massvtx);
          fgmassCut[2][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < dzbin; ++icut)
        {// dz
          double cutdz = 1*(icut+1);
          if (TMath::Abs(dzed)<cutdz) fgmassDz[icut][icent_ee]->Fill(massvtx, phoT_std);
        }

        //=====================================================================
        //         shrink the sample size and perform reconstruction
        //=====================================================================
        if ( !dz_cut || !phiv_cut ) continue;
	cout<<"creating my pair"<<endl;
        CNTDE mypair;
        reco.findIntersection(&mytrk1, &mytrk2, &mypair, zVtx); // assume zvertex=zVtx
        float radius_r = mypair.GetR();
        float phi_r    = (mypair.GetEphi()+mypair.GetPphi())/2;
        float theta_r  = (mypair.GetEtheta()+mypair.GetPtheta())/2;
        float dtheta   = mypair.GetEtheta()-mypair.GetPtheta();

        TVector3 mom1_r  = reco.findMomentum(&mytrk1, radius_r, phi_r, theta_r, zVtx);
        TVector3 mom2_r  = reco.findMomentum(&mytrk2, radius_r, phi_r, theta_r, zVtx);
        float pT1_r = mom1_r.Pt();
        float pT2_r = mom2_r.Pt();

        float r_cut = 0, phi_cut = 0, dtheta_cut = 0;
        if ( (radius_r>=9) && (radius_r<=23) ) r_cut = 1;
        // if ( ((phi_r>=-1) && (phi_r<=1)) || ((phi_r>=2) && (phi_r<=4)) ) phi_cut = 1;
        phi_cut = 1;
        if ( TMath::Abs(dtheta)<=DTHETA ) dtheta_cut = 1;

        //=====================================================
        //        implicit cut for solution of r_{conv}
        //=====================================================
        if ( radius_r<0 ) continue; // solution cut, including out of range case
        float phoT = getPT(pT1*cos(phi_r), pT1*sin(phi_r), pz1, pT2*cos(phi_r), pT2*sin(phi_r), pz2); // pz can be replaced with pT1/cot(theta_r)

        fgmassCutMBAllpt[3]->Fill(massvtx);
        fgmassCutAllpt[3][icent_ee]->Fill(massvtx);
        fgmassCut[3][icent_ee]->Fill(massvtx, phoT_std);        
        
        //=====================================================
        //                    r_{conv} cut
        //=====================================================
        fgRadius[icent_ee]->Fill(radius_r, phoT_std);
        if ( r_cut )
        {
          fgmassCutMBAllpt[4]->Fill(massvtx);
          fgmassCutAllpt[4][icent_ee]->Fill(massvtx);
          fgmassCut[4][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < rbin1; ++icut)
        {// radius1
          double cutradius1 = 0.5*(icut+1);
          if ((radius_r>=16)&&(radius_r<16+cutradius1)) fgmassRadius1[icut][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < rbin2; ++icut)
        {// radius2
          double cutradius2 = 0.5*(icut+1);
          if ((radius_r<16)&&(radius_r>=16-cutradius2)) fgmassRadius2[icut][icent_ee]->Fill(massvtx, phoT_std);
        }

        //=====================================================
        //                   phi_{conv} cut
        //=====================================================
        fgPhi[icent_ee]->Fill(phi_r, phoT_std);
        if ( phi_cut )
        {
          fgmassCutMBAllpt[5]->Fill(massvtx);
          fgmassCutAllpt[5][icent_ee]->Fill(massvtx);
          fgmassCut[5][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < phibin1; ++icut)
        {// phi1
          double cutphi1 = 0.05*(icut+1);
          if ((phi_r>=0.3)&&(phi_r<0.3+cutphi1)) fgmassPhi1[icut][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < phibin2; ++icut)
        {// phi2
          double cutphi2 = 0.05*(icut+1);
          if ((phi_r<0.3)&&(phi_r>=0.3-cutphi2)) fgmassPhi2[icut][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < phibin3; ++icut)
        {// phi3
          double cutphi3 = 0.05*(icut+1);
          if ((phi_r>=3.2)&&(phi_r<3.2+cutphi3)) fgmassPhi3[icut][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < phibin4; ++icut)
        {// phi4
          double cutphi4 = 0.05*(icut+1);
          if ((phi_r<3.2)&&(phi_r>=3.2-cutphi4)) fgmassPhi4[icut][icent_ee]->Fill(massvtx, phoT_std);
        }

        //=====================================================
        //                      dtheta cut
        //=====================================================
        fgTheta[icent_ee]->Fill(theta_r, phoT_std);
        fgDtheta[icent_ee]->Fill(dtheta, phoT_std);
        if ( dtheta_cut )
        {
          fgmassCutMBAllpt[6]->Fill(massvtx);
          fgmassCutAllpt[6][icent_ee]->Fill(massvtx);
          fgmassCut[6][icent_ee]->Fill(massvtx, phoT_std);
        }
        for(int icut = 0; icut < dthetabin; ++icut)
        {// dz
          double cutdtheta = 0.001*(icut+1);
          if (TMath::Abs(dtheta)<cutdtheta) fgmassDtheta[icut][icent_ee]->Fill(massvtx, phoT_std);
        }
        
        //=====================================================
        //                      all cuts  
        //=====================================================
        if ( !r_cut || !phi_cut || !dtheta_cut ) continue;
        fgmassCutMBAllpt[7]->Fill(massvtx);
        fgmassCutAllpt[7][icent_ee]->Fill(massvtx);
        fgmassCut[7][icent_ee]->Fill(massvtx, phoT_std);

        //====================================================================
        //                     find converted photon
        //====================================================================
        fgPhotonPt->Fill(phoT);
        fgStdMass->Fill(massvtx);
        fgCStdMass->Fill(getMass(pT1*cos(phi_r), pT1*sin(phi_r), pz1, pT2*cos(phi_r), pT2*sin(phi_r), pz2));
        fgRecoMass->Fill(getMass(pT1_r*cos(phi_r), pT1_r*sin(phi_r), pz1, pT2_r*cos(phi_r), pT2_r*sin(phi_r), pz2));
      }
    }// jump out from ireal1 loop
  }
	    
  cout<<"writing fg histograms to file"<<endl;
  TFile* f2 = new TFile(Form("mixingfg%d.root", job),"RECREATE");
 
  fgPhiDC->Write();
  
  fgPhotonPt->Write();
  fgStdMass->Write();
	fgCStdMass->Write();
	fgRecoMass->Write(); 

  for (int i = 0; i < centbin; ++i)
  {
    fgPhiv[i]->Write(Form("fgPhiv%d",i));
    fgRadius[i]->Write(Form("fgRadius%d",i));
    fgPhi[i]->Write(Form("fgPhi%d",i));
    fgTheta[i]->Write(Form("fgTheta%d",i));
    fgDz[i]->Write(Form("fgDz%d",i));
    fgDtheta[i]->Write(Form("fgDtheta%d",i));
    // fgDphi[i]->Write(Form("fgDphi%d",i));
  }
  
  for (int i = 0; i < cutbin; ++i)
  {
    fgmassCutMBAllpt[i]->Write(Form("fgmassCutMBAllpt%d", i));
    for (int j = 0; j < centbin; ++j)
    {
      fgmassCutAllpt[i][j]->Write(Form("fgmassCutAllpt%d_%d", i, j)); // 1D hist
      fgmassCut[i][j]->Write(Form("fgmassCut%d_%d", i, j)); // 2D hist
    }
  }

  for(int i = 0; i < phivbin; ++i)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassPhiv[i][j]->Write(Form("fgmassPhiv%d_%d", i, j));
    }
  }
  for(int i = 0; i < dzbin; ++i)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassDz[i][j]->Write(Form("fgmassDz%d_%d", i, j));
    }
  }
  for(int i = 0; i < rbin1; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassRadius1[i][j]->Write(Form("fgmassRadius1%d_%d", i, j));
    }
  }
  for(int i = 0; i < rbin2; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassRadius2[i][j]->Write(Form("fgmassRadius2%d_%d", i, j));
    }
  }
  for(int i = 0; i < phibin1; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassPhi1[i][j]->Write(Form("fgmassPhi1%d_%d", i, j));
    }
  }
  for(int i = 0; i < phibin2; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassPhi2[i][j]->Write(Form("fgmassPhi2%d_%d", i, j));
    }
  }
  for(int i = 0; i < phibin3; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassPhi3[i][j]->Write(Form("fgmassPhi3%d_%d", i, j));
    }
  }
  for(int i = 0; i < phibin4; i++)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassPhi4[i][j]->Write(Form("fgmassPhi4%d_%d", i, j));
    }
  }
  for(int i = 0; i < dthetabin; ++i)
  {
    for (int j = 0; j < centbin; ++j)
    {
      fgmassDtheta[i][j]->Write(Form("fgmassDtheta%d_%d", i, j));
    }
  }

  f2->Write();

  f1->Close();
  f2->Close();

  delete f1;
  delete f2;

  return 0;
}
