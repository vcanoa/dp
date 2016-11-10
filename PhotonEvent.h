#ifndef __PHOTONEVENT_H__
#define __PHOTONEVENT_H__

#include "CNTE.h"
#include "CNTDE.h"
#include "EMCC.h"

#include <vector>

using namespace std;

class PhotonEvent {
 private:
  float fBBCcharge;
  float fZVertex;
  float fCentrality;
  float fPsi_BBC;
  float fPsi_BBCN;
  float fPsi_BBCS;
  float fSvxz;
  float fRun;
  vector<CNTE> fP;
  vector<CNTE> fN;
  vector<CNTDE> fDE;
  vector<EMCC> fEMCC;

 public:
  PhotonEvent();
  PhotonEvent(const PhotonEvent &source);
  virtual ~PhotonEvent();
  
  void CloneFrom(PhotonEvent *source);
  void Clear();
  
  void SetRunNumber(int val)    { fRun = val; };
  void SetBBCcharge(float val)  { fBBCcharge = val; };
  void SetVtxZ(float val)       { fZVertex = val; }
  void SetCentrality(float val) { fCentrality = val; };
  void SetPsiBBC(float val)     { fPsi_BBC = val; };
  void SetPsiBBCN(float val)    { fPsi_BBCN = val; };
  void SetPsiBBCS(float val)    { fPsi_BBCS = val; };
  void SetSvxZ(float val)       { fSvxz = val; };

  int   GetRunNumber()  { return fRun; }
  float GetBBCcharge()  { return fBBCcharge; };
  float GetVtxZ()       { return fZVertex; }
  float GetCentrality() { return fCentrality; };
  float GetPsiBBC()     { return fPsi_BBC; };
  float GetPsiBBCN()    { return fPsi_BBCN; };
  float GetPsiBBCS()    { return fPsi_BBCS; };
  float GetSvxZ()       { return fSvxz; };
  
  void AddPTrack(CNTE val)  { fP.push_back(val); };
  void AddNTrack(CNTE val)  { fN.push_back(val); };
  void AddPair(CNTDE val)   { fDE.push_back(val); };
  void AddCluster(EMCC val) { fEMCC.push_back(val); };
  
  int GetNPtracks()  { return fP.size(); };
  int GetNEtracks()  { return fN.size(); };
  int GetNpairs()    { return fDE.size(); };
  int GetNclusters() { return fEMCC.size(); };
  
  CNTE*  GetPtrack(int i) { return &(fP[i]); };
  CNTE*  GetNtrack(int i) { return &(fN[i]); };
  EMCC*  GetCluster(int i){ return &(fEMCC[i]); };
  CNTDE* GetPair(int i)   { return &(fDE[i]); };
  
  vector<CNTE>  GetPTracks()  { return fP; };
  vector<CNTE>  GetNTracks()  { return fN; };
  vector<EMCC>  GetClusters() { return fEMCC; };
  vector<CNTDE> GetPairs()    { return fDE; };
};

#endif
