#ifndef __CNTE_H__
#define __CNTE_H__

#include "TObject.h"

class CNTE : public TObject {
 private:
  int fDCarm;
  int fDCside;
  float fPx;
  float fPy;
  float fPz;
  float fPhi;
  float fPhi0;
  float fTheta0;
  float fZed;
  float fAlpha;
  int   fCharge;
  int   fEMCid;
  float fEcore;
  float fDisp;
  int   fN0;
  float fChi2;
  float fNpe0;
  float fDep;
  float fProb;
  float fEMCdz;
  float fEMCdphi;
  float fEMCx;
  float fEMCy;
  float fEMCz;
  float fCenterPhi;
  float fCenterZ;
  float fPPC1x;
  float fPPC1y;
  float fPPC1z;
  int   fEsect;
  int   fYsect;
  int   fZsect;
  
 public:
  CNTE();
  virtual ~CNTE();
  
  int   GetDCarm()     const { return fDCarm; };
  int   GetDCside()    const { return fDCside; };
  float GetPx()        const { return fPx; };
  float GetPy()        const { return fPy; };
  float GetPz()        const { return fPz; };
  float GetPhi()       const { return fPhi; };
  float GetPhi0()      const { return fPhi0; };
  float GetTheta0()    const { return fTheta0; };
  float GetZed()       const { return fZed; };
  float GetAlpha()     const { return fAlpha; };
  int   GetCharge()    const { return fCharge; };
  int   GetEMCid()     const { return fEMCid; };
  float GetEcore()     const { return fEcore; };
  float GetDisp()      const { return fDisp; };
  int   GetN0()        const { return fN0; };
  float GetChi2()      const { return fChi2; };
  float GetNpe0()      const { return fNpe0; };
  float GetDep()       const { return fDep; };
  float GetProb()      const { return fProb; };
  float GetEMCdz()     const { return fEMCdz; };
  float GetEMCdphi()   const { return fEMCdphi; };
  float GetEMCx()      const { return fEMCx; };
  float GetEMCy()      const { return fEMCy; };
  float GetEMCz()      const { return fEMCz; };
  float GetCenterPhi() const { return fCenterPhi; };
  float GetCenterZ()   const { return fCenterZ; };
  float GetPPC1x()     const { return fPPC1x; };
  float GetPPC1y()     const { return fPPC1y; };
  float GetPPC1z()     const { return fPPC1z; };
  int   GetEsect()     const { return fEsect; };
  int   GetYsect()     const { return fYsect; };
  int   GetZsect()     const { return fZsect; };
  
  void SetDCarm(int val)        { fDCarm = val; };
  void SetDCside(int val)       { fDCside = val; };
  void SetPx(float val)         { fPx = val; };
  void SetPy(float val)         { fPy = val; };
  void SetPz(float val)         { fPz = val; };
  void SetPhi(float val)        { fPhi = val; };
  void SetZed(float val)        { fZed = val; };
  void SetPhi0(float val)       { fPhi0 = val; };
  void SetTheta0(float val)     { fTheta0 = val; };
  void SetAlpha(float val)      { fAlpha = val; };
  void SetCharge(int val)       { fCharge = val; };
  void SetEMCid(int val)        { fEMCid = val; };
  void SetEcore(float val)      { fEcore = val; };
  void SetDisp(float val)       { fDisp = val; };
  void SetN0(int val)           { fN0 = val; };
  void SetChi2(float val)       { fChi2 = val; };
  void SetNpe0(float val)       { fNpe0 = val; };
  void SetDep(float val)        { fDep = val; };
  void SetProb(float val)       { fProb = val; };
  void SetEMCdz(float val)      { fEMCdz = val; };
  void SetEMCdphi(float val)    { fEMCdphi = val; }
  void SetEMCx(float val)       { fEMCx = val; };
  void SetEMCy(float val)       { fEMCy = val; };
  void SetEMCz(float val)       { fEMCz = val; };
  void SetCenterPhi(float val)  { fCenterPhi = val; };
  void SetCenterZ(float val)    { fCenterZ = val; };
  void SetPPC1x(float val)      { fPPC1x = val; };
  void SetPPC1y(float val)      { fPPC1y = val; };
  void SetPPC1z(float val)      { fPPC1z = val; };
  void SetEsect(float val)      { fEsect = val; };
  void SetYsect(float val)      { fYsect = val; };
  void SetZsect(float val)      { fZsect = val; };
  
  ClassDef(CNTE,1)
};
#endif
