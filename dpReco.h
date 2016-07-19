#ifndef __DPRECO_H__
#define __DPRECO_H__

#include <vector>

#include "TVector3.h"
#include "TMatrixD.h"

#include "CNTE.h"
#include "CNTDE.h"

using namespace std;

//======================
class P33 {
 public:
  P33() {}
  P33(float X, float Y, float Z, float A, float B, float C);
  P33( const P33 &other );
  float fX;
  float fY;
  float fZ;
  float fDep[3];
};
//======================
class P33_S {
 public:
  P33_S() {}
  P33_S(float X_LO, float X_HI, float Y_LO, float Y_HI,
	float Z_LO, float Z_HI, int MLEV, int LEV=0);
  virtual ~P33_S();
  virtual bool Insert( P33 const& point );
  virtual void AppendList(vector<P33>& point_list,
			  float X_LO, float X_HI, float Y_LO, float Y_HI,
			  float Z_LO, float Z_HI);
  float fX_lo;
  float fX_hi;
  float fY_lo;
  float fY_hi;
  float fZ_lo;
  float fZ_hi;
  int fLevel;
  int fMaxLevel;
  P33_S *fContainers[2][2][2];
};
//======================
class P33_SE: public P33_S {
 public:
  P33_SE() {}
  P33_SE( float X_LO, float X_HI, float Y_LO, float Y_HI,
	  float Z_LO, float Z_HI, int MLEV, int LEV=0 );
  virtual ~P33_SE();
  virtual bool Insert( P33 const& point );
  virtual void AppendList( vector<P33>& point_list,
			   float X_LO, float X_HI, float Y_LO, float Y_HI,
			   float Z_LO, float Z_HI );
  vector<P33> fPoints;
};
//======================
class LUT {
 public:
  LUT() {}
  LUT( const char* lookup_file_name );
  P33_S alpha_r_z;
};
//======================
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
		      vector<P33>& points_temp,
		      vector<P33>& points );
  TMatrixD PerformFit(int ind, vector<P33> const& points );
  float LookupFit(float alpha, float r, float z_DC_m_z_ver, int ind, 
		  P33_S& alpha_r_z, float& lookup_ind );
  bool ProjectPhi(float alpha, float phiDC,
		  float r, float z_DC_m_z_ver,
		  P33_S& alpha_r_z, float& phi_r );
  bool ProjectTheta(float alpha, float r, float z_DC_m_z_ver,
		    P33_S& alpha_r_z, float& theta_r );
  bool DeltaPhi(float alpha_e, float alpha_p,
		float phi_e, float phi_p,
		float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
		P33_S& alpha_r_z, float& dphi );
  bool FindR(float alpha_e, float alpha_p,
	     float phi_e, float phi_p,
	     float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
	     P33_S& alpha_r_z, float& r_conv );

  LUT* fLUT;
};

#endif /*__DPRECO_H__*/

