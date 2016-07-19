#ifndef __PHPhotonEVENT_H__
#define __PHPhotonEVENT_H__

#include <vector>
#include "PHObject.h"
#include "PhotonEvent.h"

class PHPhotonEvent:public PHObject {
 public:
  PHPhotonEvent() {}
  virtual ~PHPhotonEvent() {fEvent.Clear();}
  PhotonEvent* GetEvent() {return &fEvent;}
 private:
  PhotonEvent fEvent;
};
#endif /* __PHPhotonEVENT_H__ */
