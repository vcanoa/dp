#ifndef __CNTDE_H__
#define __CNTDE_H__

class CNTDE {
 private:
  float fEphi;
  float fPphi;
  float fEtheta;
  float fPtheta;
  float fR;
  int fEid;
  int fPid;
  
 public:
  CNTDE();
  virtual ~CNTDE();
  
  float GetEphi()   const { return fEphi; }
  float GetPphi()   const { return fPphi; }
  float GetEtheta() const { return fEtheta; }
  float GetPtheta() const { return fPtheta; }
  float GetR()      const { return fR; }
  float GetEid()    const { return fEid; }
  float GetPid()    const { return fPid; }
  
  void SetEphi(float val)   { fEphi = val; }
  void SetPphi(float val)   { fPphi = val; }
  void SetEtheta(float val) { fEtheta = val; }
  void SetPtheta(float val) { fPtheta = val; }
  void SetR(float val)      { fR = val; }
  void SetEid(int val)      { fEid = val; }
  void SetPid(int val )     { fPid = val; }
  
};

#endif
