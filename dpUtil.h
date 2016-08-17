#ifndef __DPUTIL_H__
#define __DPUTIL_H__

#include <vector>

#include "TVector3.h"
#include "TMatrixD.h"

using namespace std;

//======================
class dpOBJA {
 public:
  dpOBJA() {}
  dpOBJA(float X, float Y, float Z, float A, float B, float C);
  dpOBJA( const dpOBJA &other );
  float fX;
  float fY;
  float fZ;
  float fDep[3];
};
//======================
class dpOBJB {
 public:
  dpOBJB() {}
  dpOBJB(float X_LO, float X_HI, float Y_LO, float Y_HI,
	 float Z_LO, float Z_HI, int MLEV, int LEV=0);
  virtual ~dpOBJB();
  virtual bool Insert( dpOBJA const& point );
  virtual void AppendList(vector<dpOBJA>& point_list,
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
  dpOBJB *fContainers[2][2][2];
};
//======================
class dpOBJAB: public dpOBJB {
 public:
  dpOBJAB() {}
  virtual ~dpOBJAB() {}
  dpOBJAB( float X_LO, float X_HI, float Y_LO, float Y_HI,
	   float Z_LO, float Z_HI, int MLEV, int LEV=0 );
  virtual bool Insert( dpOBJA const& point );
  virtual void AppendList( vector<dpOBJA>& point_list,
			   float X_LO, float X_HI, float Y_LO, float Y_HI,
			   float Z_LO, float Z_HI );
  vector<dpOBJA> fPoints;
};
//======================
class LUT {
 public:
  LUT() {}
  LUT( const char* lookup_file_name );
  dpOBJB alpha_r_z;
};
#endif /*__DPUTIL_H__*/

