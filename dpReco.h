#ifndef __DPRECO_H__
#define __DPRECO_H__

#include <vector>

#include "TVector3.h"
#include "TMatrixD.h"

#include "CNTE.h"
#include "CNTDE.h"
#include "dpUtil.h"

using namespace std;

class dpReco {
 public:
  float fMaxR;
  float fLookupAlphaDeltaBase;
  float fLookupRDeltaBase;
  float fLookupZDeltaBase;
  dpReco() {}
  dpReco(dpReco const&);
  dpReco( const char* lookup_file_name );
  ~dpReco();

  void findIntersection(CNTE *trk1, CNTE *trk2, CNTDE *pair, float zvertex);
  TVector3 findMomentum(CNTE *trk, float r_conv, float phi_conv,
			float theta_conv, float zvertex);
 private:
  float EvaluateFit(TMatrixD& beta, float alpha, float r, float z_DC_m_z_ver );
  void RemoveOutliers(int ind, TMatrixD& beta,
		      vector<dpOBJA>& points_temp,
		      vector<dpOBJA>& points );
  TMatrixD PerformFit(int ind, vector<dpOBJA> const& points );
  float LookupFit(float alpha, float r, float z_DC_m_z_ver, int ind, 
		  dpOBJB& alpha_r_z, float& lookup_ind );
  bool ProjectPhi(float alpha, float phiDC,
		  float r, float z_DC_m_z_ver,
		  dpOBJB& alpha_r_z, float& phi_r );
  bool ProjectTheta(float alpha, float r, float z_DC_m_z_ver,
		    dpOBJB& alpha_r_z, float& theta_r );
  bool DeltaPhi(float alpha_e, float alpha_p,
		float phi_e, float phi_p,
		float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
		dpOBJB& alpha_r_z, float& dphi );
  bool FindR(float alpha_e, float alpha_p,
	     float phi_e, float phi_p,
	     float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
	     dpOBJB& alpha_r_z, float& r_conv );

  LUT* fLUT;
};

#endif /*__DPRECO_H__*/

