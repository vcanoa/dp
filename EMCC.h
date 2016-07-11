#ifndef __EMCC_H__
#define __EMCC_H__

#include "TObject.h"

class EMCC : public TObject {
 private:
  int   fArm;
  int   fID;
  float fX;
  float fY;           
  float fZ;           
  float fEcore;
  float fProbPhoton;
  float fEMCdz;
  float fEMCdphi;
  
 public:
  EMCC();
  virtual ~EMCC();
  
  int      GetArm() const { return fArm; };
  int      GetID() const { return fID; };
  float    GetX() const { return fX; };
  float    GetY() const { return fY; };
  float    GetZ() const { return fZ; };
  float    GetEcore() const { return fEcore; };
  float    GetProbPhoton() const { return fProbPhoton; };
  float    GetEMCdz() const { return fEMCdz; };
  float    GetEMCdphi() const { return fEMCdphi; };
  
  
  void     SetArm(int val) { fArm = val; }
  void     SetID(int val) { fID = val; }
  void     SetX(float val) { fX = val; };
  void     SetY(float val) { fY = val; };
  void     SetZ(float val) { fZ = val; };
  void     SetEcore(float val) { fEcore = val; };
  void     SetProb(float val) { fProbPhoton = val; };
  void     SetEmcdz(float val) { fEMCdz = val; };
  void     SetEmcdphi(float val) { fEMCdphi = val; };
  
  ClassDef(EMCC,1)
};

#endif
