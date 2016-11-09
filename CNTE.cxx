#include "CNTE.h"

CNTE::CNTE() {
}
CNTE::CNTE(const CNTE &source){
  fDCarm=source.GetDCarm();
  fDCside=source.GetDCside();
  fPx=source.GetPx();
  fPy=source.GetPy();
  fPz=source.GetPz();
  fPhi=source.GetPhi();
  fPhi0=source.GetPhi0();
  fTheta0=source.GetTheta0();
  fZed=source.GetZed();
  fAlpha=source.GetAlpha();
  fCharge=source.GetCharge();
  fEMCid=source.GetEMCid();
  fEcore=source.GetEcore();
  fDisp=source.GetDisp();
  fN0=source.GetN0();
  fChi2=source.GetChi2();
  fNpe0=source.GetNpe0();
  fDep=source.GetDep();
  fProb=source.GetProb();
  fEMCdz=source.GetEMCdz();
  fEMCdphi=source.GetEMCdphi();
  fEMCx=source.GetEMCx();
  fEMCy=source.GetEMCy();
  fEMCz=source.GetEMCz();
  fCenterPhi=source.GetCenterPhi();
  fCenterZ=source.GetCenterZ();
  fPPC1x=source.GetPPC1x();
  fPPC1y=source.GetPPC1y();
  fPPC1z=source.GetPPC1z();
  fEsect=source.GetEsect();
  fYsect=source.GetYsect();
  fZsect=source.GetZsect();
}
CNTE::CNTE(CNTE *source){
  fDCarm=source->GetDCarm();
  fDCside=source->GetDCside();
  fPx=source->GetPx();
  fPy=source->GetPy();
  fPz=source->GetPz();
  fPhi=source->GetPhi();
  fPhi0=source->GetPhi0();
  fTheta0=source->GetTheta0();
  fZed=source->GetZed();
  fAlpha=source->GetAlpha();
  fCharge=source->GetCharge();
  fEMCid=source->GetEMCid();
  fEcore=source->GetEcore();
  fDisp=source->GetDisp();
  fN0=source->GetN0();
  fChi2=source->GetChi2();
  fNpe0=source->GetNpe0();
  fDep=source->GetDep();
  fProb=source->GetProb();
  fEMCdz=source->GetEMCdz();
  fEMCdphi=source->GetEMCdphi();
  fEMCx=source->GetEMCx();
  fEMCy=source->GetEMCy();
  fEMCz=source->GetEMCz();
  fCenterPhi=source->GetCenterPhi();
  fCenterZ=source->GetCenterZ();
  fPPC1x=source->GetPPC1x();
  fPPC1y=source->GetPPC1y();
  fPPC1z=source->GetPPC1z();
  fEsect=source->GetEsect();
  fYsect=source->GetYsect();
  fZsect=source->GetZsect();
}
CNTE::~CNTE() {
}
