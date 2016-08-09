#include "dpUtil.h"
#include <TFile.h>
#include <TNtuple.h>
#include <TMatrixD.h>
#include <iostream>
#include <TMath.h>

using namespace std;

//======================
dpOBJA::dpOBJA(float X, float Y, float Z, float A, float B, float C) : fX(X), fY(Y), fZ(Z) {
  fDep[0] = A;
  fDep[1] = B;
  fDep[2] = C;
}
dpOBJA::dpOBJA(const dpOBJA &other ) {
  fX = other.fX;
  fY = other.fY;
  fZ = other.fZ;
  fDep[0] = other.fDep[0];
  fDep[1] = other.fDep[1];
  fDep[2] = other.fDep[2];
}
//======================
dpOBJB::dpOBJB(float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV ) : 
  fX_lo(X_LO), fX_hi(X_HI), fY_lo(Y_LO), fY_hi(Y_HI), fZ_lo(Z_LO), fZ_hi(Z_HI), fLevel(LEV), fMaxLevel(MLEV) {
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)
      for(unsigned int k=0;k<2;++k)
	fContainers[i][j][k]=NULL;
}
dpOBJB::~dpOBJB() {
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)
      for(unsigned int k=0;k<2;++k)
	if(fContainers[i][j][k]!=NULL)
	  delete fContainers[i][j][k];
}
void dpOBJB::AppendList(vector<dpOBJA>& point_list, float X_LO, float X_HI,
		       float Y_LO, float Y_HI, float Z_LO, float Z_HI ) {
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)
      for(unsigned int k=0;k<2;++k) {
	if(fContainers[i][j][k]==NULL)
	  continue;
	if( (fContainers[i][j][k]->fX_hi<X_LO) ||
	    (fContainers[i][j][k]->fX_lo>X_HI) ||
	    (fContainers[i][j][k]->fY_hi<Y_LO) ||
	    (fContainers[i][j][k]->fY_lo>Y_HI) ||
	    (fContainers[i][j][k]->fZ_hi<Z_LO) ||
	    (fContainers[i][j][k]->fZ_lo>Z_HI) )
	  continue;
	fContainers[i][j][k]->AppendList( point_list, X_LO, X_HI, Y_LO, Y_HI, Z_LO, Z_HI );
      }
}
//======================
dpOBJAB::dpOBJAB( float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV ) : 
  dpOBJB( X_LO,X_HI,Y_LO,Y_HI,Z_LO,Z_HI,MLEV,LEV ) {
}

bool dpOBJAB::Insert( dpOBJA const& point ) {
  fPoints.push_back(point);
  return true;
}
void dpOBJAB::AppendList(vector<dpOBJA>& point_list,float X_LO,float X_HI,
			float Y_LO,float Y_HI,float Z_LO,float Z_HI) {
  for(unsigned int i=0;i<fPoints.size();++i) {
    if( (fPoints[i].fX<X_LO) ||
	(fPoints[i].fX>X_HI) ||
	(fPoints[i].fY<Y_LO) ||
	(fPoints[i].fY>Y_HI) ||
	(fPoints[i].fZ<Z_LO) ||
	(fPoints[i].fZ>Z_HI) ) 
      continue;
    point_list.push_back( fPoints[i] );
  }
}
//======================
bool dpOBJB::Insert(dpOBJA const& point) {
  if( (point.fX < fX_lo) ||
      (point.fY < fY_lo) ||
      (point.fZ < fZ_lo) ||
      (point.fX > fX_hi) ||
      (point.fY > fY_hi) ||
      (point.fZ > fZ_hi) )
    return false;
  int x_ind = 0;
  int y_ind = 0;
  int z_ind = 0;
  if(point.fX > (fX_lo + 0.5*(fX_hi-fX_lo)))
    x_ind=1;
  if(point.fY > (fY_lo + 0.5*(fY_hi-fY_lo)))
    y_ind=1;
  if(point.fZ > (fZ_lo + 0.5*(fZ_hi-fZ_lo)))
    z_ind=1;
  if( fContainers[x_ind][y_ind][z_ind] == NULL ) {
    float x_lo_new = fX_lo + (float(x_ind))*0.5*(fX_hi-fX_lo);
    float x_hi_new = x_lo_new + 0.5*(fX_hi-fX_lo);
    float y_lo_new = fY_lo + (float(y_ind))*0.5*(fY_hi-fY_lo);
    float y_hi_new = y_lo_new + 0.5*(fY_hi-fY_lo);
    float z_lo_new = fZ_lo + (float(z_ind))*0.5*(fZ_hi-fZ_lo);
    float z_hi_new = z_lo_new + 0.5*(fZ_hi-fZ_lo);
    if(fLevel < fMaxLevel)
      fContainers[x_ind][y_ind][z_ind] = new dpOBJB( x_lo_new,
						    x_hi_new,
						    y_lo_new,
						    y_hi_new,
						    z_lo_new,
						    z_hi_new,
						    fMaxLevel,
						    fLevel+1 );
    else
      fContainers[x_ind][y_ind][z_ind] = new dpOBJAB(x_lo_new,
						    x_hi_new,
						    y_lo_new,
						    y_hi_new,
						    z_lo_new,
						    z_hi_new,
						    fMaxLevel,
						    fLevel+1 );
  }
  return fContainers[x_ind][y_ind][z_ind]->Insert( point );
}
//======================
LUT::LUT( const char* lookup_file_name ) : alpha_r_z( 0., 0.7, 0., 150., 0., 100., 10 ) {
  TFile f(lookup_file_name);
  TNtuple* t=0;
  f.GetObject("alpha_p_r_phi_theta_z", t);
  float alpha, p, r, phi, theta, z;
  t->SetBranchAddress("alpha",&alpha);
  t->SetBranchAddress("p",&p);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("phi",&phi);
  t->SetBranchAddress("theta",&theta);
  t->SetBranchAddress("z",&z);
  for(int i=0;i<t->GetEntries();i+=1) {
    t->GetEntry(i);
    if ( z < 0.) continue;
    if ( alpha > 0. ) continue;
    if ( phi < 0. || phi > TMath::Pi() ) continue;
    if ( theta < 0. ) continue;
    dpOBJA point(TMath::Abs(alpha),r,TMath::Abs(z), TMath::Abs(phi), p, TMath::Abs(theta));
    alpha_r_z.Insert(point);
  }
  f.Close();
}
