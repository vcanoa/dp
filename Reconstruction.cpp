#include "Reconstruction.h"
#include <TFile.h>
#include <TNtuple.h>
#include <TMatrixD.h>
#include <iostream>
#include <TMath.h>

using namespace std;

//======================
P33::P33(float X, float Y, float Z, float A, float B, float C) : fX(X), fY(Y), fZ(Z) {
  fDep[0] = A;
  fDep[1] = B;
  fDep[2] = C;
}
P33::P33(P33 const& other ) {
  fX = other.fX;
  fY = other.fY;
  fZ = other.fZ;
  fDep[0] = other.fDep[0];
  fDep[1] = other.fDep[1];
  fDep[2] = other.fDep[2];
}
//======================
P33_S::P33_S(float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV=0 ) : 
  fX_lo(X_LO), fX_hi(X_HI), fY_lo(Y_LO), fY_hi(Y_HI), fZ_lo(Z_LO), fZ_hi(Z_HI), fLevel(LEV), fMaxLevel(MLEV) {
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)
      for(unsigned int k=0;k<2;++k)
	fContainers[i][j][k]=NULL;
}
//======================
P33_S::~P33_S() {
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)
      for(unsigned int k=0;k<2;++k)
	if(fContainers[i][j][k]!=NULL)
	  delete fContainers[i][j][k];
}
//======================
void P33_S::AppendList( vector<PointValue_3x3>& point_list, float X_LO, float X_HI,
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

bool P33_S::Insert( PointValue_3x3 const& point ) {
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
    if(level < maxlevel)
      fContainers[x_ind][y_ind][z_ind] = new P33_S( x_lo_new,
						    x_hi_new,
						    y_lo_new,
						    y_hi_new,
						    z_lo_new,
						    z_hi_new,
						    maxlevel,
						    level+1 );
    else
      fContainers[x_ind][y_ind][z_ind] = new P33_SE(x_lo_new,
						    x_hi_new,
						    y_lo_new,
						    y_hi_new,
						    z_lo_new,
						    z_hi_new,
						    maxlevel,
						    level+1 );
  }
  return fContainers[x_ind][y_ind][z_ind]->Insert( point );
}

bool P33_SE:Insert( P33 const& point ) {
  fPoints.push_back(point);
  return true;
}

void P33End::AppendList(vector<P33>& point_list,float X_LO,float X_HI,
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
    PointValue_3x3 point(TMath::Abs(alpha),r,TMath::Abs(z), TMath::Abs(phi), p, TMath::Abs(theta));
    alpha_r_z.Insert(point);
  }
  f.Close();
}

Reconstruction::Reconstruction( const char* lookup_file_name ) {
  fLUT = new LUT(lookup_file_name);
}

Reconstruction::~Reconstruction() {
  delete fLUT;
}

float Reconstruction::EvaluateFit( TMatrixD& beta, float alpha, float r, float z_DC_m_z_ver ) {
  float retval = 0.;
  retval += beta(0,0);
  retval += beta(1,0)*alpha;
  retval += beta(2,0)*r;
  retval += beta(3,0)*z_DC_m_z_ver;
  retval += beta(4,0)*alpha*alpha;
  retval += beta(5,0)*alpha*r;
  retval += beta(6,0)*alpha*z_DC_m_z_ver;
  retval += beta(7,0)*r*r;
  retval += beta(8,0)*r*z_DC_m_z_ver;
  retval += beta(9,0)*z_DC_m_z_ver*z_DC_m_z_ver;
  
  return retval;<
}

void Reconstruction::RemoveOutliers( int ind, TMatrixD& beta, vector<PointValue_3x3>& points_temp, vector<PointValue_3x3>& points ) {
  points.clear();

  float sig = 0.;
  for( unsigned int i=0; i<points_temp.size(); ++i )
    {
      float fval = EvaluateFit( beta, points_temp[i].x, points_temp[i].y, points_temp[i].z );
      sig += ( fval - points_temp[i].dep[ind] )*( fval - points_temp[i].dep[ind] );
    }
  sig /= ((float)(points_temp.size()));
  sig = sqrt(sig);
  
  for( unsigned int i=0; i<points_temp.size(); ++i )
    {
      float fval = EvaluateFit( beta, points_temp[i].x, points_temp[i].y, points_temp[i].z );
      if( fabs( fval - points_temp[i].dep[ind] ) < 3.*sig ){points.push_back(points_temp[i]);}
    }
}

TMatrixD Reconstruction::PerformFit( int ind, vector<PointValue_3x3> const& points )
{
  TMatrixD X( points.size(), 10 );
  TMatrixD y( points.size(), 1  );
  for(unsigned int i=0;i<points.size();i+=1)
    {
      y(i,0) = points[i].dep[ind];
      X(i,0) = 1.;
      X(i,1) = points[i].x;
      X(i,2) = points[i].y;
      X(i,3) = points[i].z;
      X(i,4) = points[i].x * points[i].x;
      X(i,5) = points[i].x * points[i].y;
      X(i,6) = points[i].x * points[i].z;
      X(i,7) = points[i].y * points[i].y;
      X(i,8) = points[i].y * points[i].z;
      X(i,9) = points[i].z * points[i].z;
    }
  
  TMatrixD Xt( 10, points.size() );
  Xt.Transpose(X);
  
  TMatrixD XtX = Xt * X;
  XtX.Invert();
  
  TMatrixD beta = XtX * (Xt * y);
  
  return beta;
}


// 0 <--> phi
// 1 <--> p
// 2 <--> theta
float Reconstruction::LookupFit( float alpha, float r, float z_DC_m_z_ver, int ind, P33& alpha_r_z, float& lookup_ind ) {
  // r can be -9999 when find_conv() fails
  if ( r<0 || r>max_r ) { lookup_ind = -9999.; return 0; }
  if ( alpha<0 || z_DC_m_z_ver<0 ) {
    cout<<"input error!!!"<<endl;
    cout<<"lookup table contain only e- entries, please convert e+ to e- before fitting with lookup table"<<endl;
    lookup_ind = -9999;
    return 0;
  }
  float lookup_alpha_delta = lookup_alpha_delta_base;
  float lookup_r_delta = lookup_r_delta_base;
  float lookup_z_delta = lookup_z_delta_base;
  vector<PointValue_3x3> points;
  alpha_r_z.AppendList( points, alpha - lookup_alpha_delta, alpha + lookup_alpha_delta,
			 r - lookup_r_delta, r + lookup_r_delta,
			 z_DC_m_z_ver - lookup_z_delta, z_DC_m_z_ver + lookup_z_delta );
  if( points.size() < 20 ) {
    lookup_ind = -9999.;
    return false;
  }
	// p0 + p1*x + p2*y + p3*z + p4*x^2 + p5*xy + p6*xz + p7*y^2 + p8*yz + p9*z^2
  TMatrixD beta = PerformFit( ind, points );
  {
    vector<PointValue_3x3> points_temp = points;
    RemoveOutliers( ind, beta, points_temp, points );
  }
  if( points.size() < 20 ) {
    lookup_ind = -9999.;
    return false;
  }
  TMatrixD beta2 = PerformFit( ind, points );
  lookup_ind = EvaluateFit( beta2, alpha, r, z_DC_m_z_ver );
  return true;
}


bool Reconstruction::ProjectPhi( float alpha, float phiDC, float r, float z_DC_m_z_ver, P33& alpha_r_z, float& phi_r ) {
  //DC coord sys: phi -0.5Pi~1.5Pi
  if ( (r<0) || (r>max_r) ) { phi_r=-9999.; return 0; }
  float lookup_phi = 0;
  if ( !LookupFit(TMath::Abs(alpha), r, TMath::Abs(z_DC_m_z_ver), 0, alpha_r_z, lookup_phi) ) { phi_r=-9999.; return 0; }
  if ( alpha>0 ) {
    // electron!
    if ( lookup_phi<0 ) {
      cout<<"error!!! phi_conv-phi_DC can't be negative for electron using ++ field"<<endl;
      cout<<"lookup_phi "<<lookup_phi<<endl;
      cout<<"alpha "<<alpha<<endl;
      cout<<"r "<<r<<endl;
      cout<<"z_DC_m_z_ver "<<TMath::Abs(z_DC_m_z_ver)<<endl;
      phi_r = -9999.;
      return 0;
    } 	
  } else {
    // positron!
    lookup_phi = - lookup_phi;
    if ( lookup_phi>0 ) {
      cout<<"error!!! phi_conv-phi_DC can't be positive for positron using ++ field"<<endl;
      cout<<"lookup_phi "<<lookup_phi<<endl;
      cout<<"alpha "<<alpha<<endl;
      cout<<"r "<<r<<endl;
      cout<<"z_DC_m_z_ver "<<TMath::Abs(z_DC_m_z_ver)<<endl;
      phi_r = -9999.;
      return 0;
    } 
  }
  phi_r = phiDC + lookup_phi;
  // periodic boundary, to ensure phi_r from [-0.5pi, 1.5pi]
  if (phi_r>1.5*TMath::Pi())  phi_r = phi_r-2*TMath::Pi();
  if (phi_r<-0.5*TMath::Pi()) phi_r = phi_r+2*TMath::Pi();
  // else already in range, phi_r = phi_r
  return 1;
}

bool Reconstruction::ProjectTheta( float alpha, float r, float z_DC_m_z_ver,
				    P33& alpha_r_z, float& theta_r ) {
  //DC coord sys: phi -0.5Pi~1.5Pi
  if ( (r<0) || (r>max_r) ) {
    theta_r=-9999.;
    return 0;
  }
  float lookup_theta = 0;
  if ( !LookupFit(TMath::Abs(alpha), r, TMath::Abs(z_DC_m_z_ver), 2, alpha_r_z, lookup_theta) ) {
    theta_r = -9999.;
    return 0;
  }
  // cout<<"TMath::Abs(z_DC_m_z_ver) "<<TMath::Abs(z_DC_m_z_ver)<<" lookup_theta "<<lookup_theta<<endl;
  if ( lookup_theta<0 ) {
    cout<<"error!!! there should be no negative theta_conv from lookup table"<<endl;
    theta_r = -9999.;
    return 0;
  }
  if ( z_DC_m_z_ver>0 ) theta_r = lookup_theta; // top half
  if ( z_DC_m_z_ver<0 ) theta_r = TMath::Pi() - lookup_theta; // bottom half
  return 1;
}

bool Reconstruction::DeltaPhi( float alpha_e, float alpha_p,
				float phi_e, float phi_p,
				float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
				P33& alpha_r_z, float& dphi ) {
  float phir_e, phir_p;
  if( !ProjectPhi(alpha_e, phi_e, r, z_DC_m_z_ver_e, alpha_r_z, phir_e) ||
      !ProjectPhi(alpha_p, phi_p, r, z_DC_m_z_ver_p, alpha_r_z, phir_p) ) {
    dphi = -9999.;
    return 0;
  }
  dphi = phir_p - phir_e;
  // the dphi we use here should alway be |dphi|<=pi
  if ( dphi>TMath::Pi() )  dphi = 2*TMath::Pi()-dphi; // 1.3pi->0.7pi
  if ( dphi<-TMath::Pi() ) dphi = -2*TMath::Pi()-dphi; // -1.3pi->-0.7pi
  return 1;
}

bool Reconstruction::FindR( float alpha_e, float alpha_p,
			    float phi_e, float phi_p,
			    float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
			    P33& alpha_r_z, float& r_conv ) {
  // assumption: DeltaPhi(r) should be monotomic increase within [0, max_r]
  // hence |DeltaPhi(r)| should be unimodal
  // use godel section search to find the minimum of |DeltaPhi(r)| within range [0, max_r]
  double a = 0, b = max_r;
  double k = (sqrt(5.) - 1.) / 2.;
  double xL = b - k * (b - a);
  double xR = a + k * (b - a);
  double EPSILON = 0.001;
  int iter = 0, max_iter = 50;
  while ( fabs(b-a) > EPSILON*(fabs(xR)+fabs(xL)) ) {
    if ( iter>max_iter ) {
      cout<<"nothing found within interation!!!"<<endl;
      r_conv = -9999.;
      return 0;
    }
    float dphi_xR, dphi_xL;
    if ( !DeltaPhi( alpha_e, alpha_p, phi_e, phi_p, xR, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi_xR ) ) {
      // if DeltaPhi in xR, [xR, b] will be the new searching range
      a = xR;
      xL = b - k * (b - a);
      xR = a + k * (b - a);
      ++iter;
      continue;
    }
    if ( !DeltaPhi( alpha_e, alpha_p, phi_e, phi_p, xL, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi_xL ) ) {
      // if DeltaPhi in xL, [xL, b] will be the new searching range, xR will be the new xL
      a = xL;
      xL = xR;
      xR = a + k * (b - a);
      ++iter;
      continue;
    }
    if( TMath::Abs(dphi_xL) < TMath::Abs(dphi_xR) ) {
      // find minimum of |dphi|
      b = xR;
      xR = xL;
      xL = b - k*(b - a);
      ++iter;
    } else {
      a = xL;
      xL = xR;
      xR = a + k * (b - a);
      ++iter;
    }
  }
  r_conv = (a + b) / 2.;
  float dphi;
  DeltaPhi( alpha_e, alpha_p, phi_e, phi_p, r_conv, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi );
  // cout<<"r_conv "<<r_conv<<" DeltaPhi "<<dphi<<endl;
  if ( dphi < -9000 ) {
    r_conv = -9999.;
    return 0;
  }
  // can't do LookupFit for the whole r range [0, max_r], dphi always -9999., and that's the minimum value
  return 1;
}

void Reconstruction::findIntersection(CNTE *pos,
				      CNTE *neg,
				      CNTDE *de,
				      float zvertex) {
  if( pos->GetAlpha()*neg->GetAlpha()>0 )
    return;

  float r_conv, phi_conv_e, phi_conv_p, theta_conv_e, theta_conv_p;
  FindR( pos->GetAlpha(), neg->GetAlpha(),
	      pos->GetPhiDC(), neg->GetPhiDC(),
	      pos->GetZDC() - zvertex, neg->GetZDC() - zvertex,
	      pimpl->alpha_r_z, r_conv );
  ProjectPhi( pos->GetAlpha(), pos->GetPhiDC(), r_conv, pos->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_e);
  ProjectPhi( neg->GetAlpha(), neg->GetPhiDC(), r_conv, neg->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_p );
  ProjectTheta( pos->GetAlpha(), r_conv, pos->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_e );
  ProjectTheta( neg->GetAlpha(), r_conv, neg->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_p );

  de->SetR( r_conv );
  if ( r_conv<-9000. ) {
    de->SetEPhi( -9999. );
    de->SetPPhi( -9999. );
    de->SetETheta( -9999. );
    de->SetPTheta( -9999. );
    return;
  }
  de->SetEphi( phi_conv_e );
  de->SetPphi( phi_conv_p );
  de->SetEtheta( theta_conv_e );
  de->SetPtheta( theta_conv_p );
}

TVector3 Reconstruction::findMomentum(CNTE* trk, float r, float phi_conv, float theta_conv, float zvertex) {
  TVector3 vec;
  float mom, phi, theta;
  LookupFit( TMath::Abs(trk->GetAlpha()), r, TMath::Abs(trk->GetZed() - zvertex), 1, pimpl->alpha_r_z, mom );
  ProjectPhi( trk->GetAlpha(), trk->GetPhi(), r, TMath::Abs(trk->GetZed() - zvertex), pimpl->alpha_r_z, phi );
  ProjectTheta( trk->GetAlpha(), r, trk->GetZDC() - zvertex, pimpl->alpha_r_z, theta );
  float pt = mom*sin(theta);
  vec.SetPtThetaPhi( pt, theta, phi );
  return vec;
}
