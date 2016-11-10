#include "PhotonEvent.h"
#include "CNTE.h"
#include "CNTDE.h"
#include "EMCC.h"

#include <vector>

using namespace std;

PhotonEvent::PhotonEvent() {
  fBBCcharge  = -999;
  fZVertex    = -999;
  fCentrality = -999;
  fPsi_BBC    = -999;
  fPsi_BBCN   = -999;
  fPsi_BBCS   = -999;
  fSvxz       = -999;
  fRun        = -999;
}

PhotonEvent::PhotonEvent(const PhotonEvent &source) {
  fBBCcharge  = source.fBBCcharge;
  fZVertex    = source.fZVertex;
  fCentrality = source.fCentrality;
  fPsi_BBC    = source.fPsi_BBC;
  fPsi_BBCN   = source.fPsi_BBCN;
  fPsi_BBCS   = source.fPsi_BBCS;
  fSvxz       = source.fSvxz;
  fRun        = source.fRun;
  for(int i=0; i!=source.fP.size(); ++i) fP.push_back( source.fP[i] );
  for(int i=0; i!=source.fN.size(); ++i) fN.push_back( source.fN[i] );
  for(int i=0; i!=source.fDE.size(); ++i) fDE.push_back( source.fDE[i] );
  for(int i=0; i!=source.fEMCC.size(); ++i) fEMCC.push_back( source.fEMCC[i] );
}

void PhotonEvent::CloneFrom(PhotonEvent *source) {
  fBBCcharge  = source->fBBCcharge;
  fZVertex    = source->fZVertex;
  fCentrality = source->fCentrality;
  fPsi_BBC    = source->fPsi_BBC;
  fPsi_BBCN   = source->fPsi_BBCN;
  fPsi_BBCS   = source->fPsi_BBCS;
  fSvxz       = source->fSvxz;
  fRun        = source->fRun;
  for(int i=0; i!=source->fP.size(); ++i) fP.push_back( source->fP[i] );
  for(int i=0; i!=source->fN.size(); ++i) fN.push_back( source->fN[i] );
  for(int i=0; i!=source->fDE.size(); ++i) fDE.push_back( source->fDE[i] );
  for(int i=0; i!=source->fEMCC.size(); ++i) fEMCC.push_back( source->fEMCC[i] );
}

PhotonEvent::~PhotonEvent() {
}

void PhotonEvent::Clear() {
  fP.clear();
  fN.clear();
  fDE.clear();
  fEMCC.clear();
}
