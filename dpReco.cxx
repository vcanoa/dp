#include "dpReco.h"
#include "dpUtil.h"
#include <TFile.h>
#include <TNtuple.h>
#include <TMatrixD.h>
#include <iostream>
#include <TMath.h>
#include "CNTE.h"
#include "CNTDE.h"

using namespace std;

dpReco::dpReco( const char* lookup_file_name ) {
  fLUT = new LUT(lookup_file_name);
  fMaxR = 30.;
  fLookupAlphaDeltaBase = 0.01;
  fLookupRDeltaBase = 2.0;
  fLookupZDeltaBase = 4.0;
}
dpReco::~dpReco() {
  delete fLUT;
}

float dpReco::EvaluateFit( TMatrixD& beta, float alpha, float r, float z_DC_m_z_ver ) {
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
  
  return retval;
}

void dpReco::RemoveOutliers( int ind, TMatrixD& beta, vector<dpOBJA>& points_temp, vector<dpOBJA>& points ) {
  points.clear();

  float sig = 0.;
  for( unsigned int i=0; i<points_temp.size(); ++i )
    {
      float fval = EvaluateFit( beta, points_temp[i].fX, points_temp[i].fY, points_temp[i].fZ );
      sig += ( fval - points_temp[i].fDep[ind] )*( fval - points_temp[i].fDep[ind] );
    }
  sig /= ((float)(points_temp.size()));
  sig = sqrt(sig);
  
  for( unsigned int i=0; i<points_temp.size(); ++i )
    {
      float fval = EvaluateFit( beta, points_temp[i].fX, points_temp[i].fY, points_temp[i].fZ );
      if( fabs( fval - points_temp[i].fDep[ind] ) < 3.*sig ){points.push_back(points_temp[i]);}
    }
}

TMatrixD dpReco::PerformFit( int ind, vector<dpOBJA> const& points )
{
  TMatrixD X( points.size(), 10 );
  TMatrixD y( points.size(), 1  );
  for(unsigned int i=0;i<points.size();i+=1)
    {
      y(i,0) = points[i].fDep[ind];
      X(i,0) = 1.;
      X(i,1) = points[i].fX;
      X(i,2) = points[i].fY;
      X(i,3) = points[i].fZ;
      X(i,4) = points[i].fX * points[i].fX;
      X(i,5) = points[i].fX * points[i].fY;
      X(i,6) = points[i].fX * points[i].fZ;
      X(i,7) = points[i].fY * points[i].fY;
      X(i,8) = points[i].fY * points[i].fZ;
      X(i,9) = points[i].fZ * points[i].fZ;
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
float dpReco::LookupFit( float alpha, float r, float z_DC_m_z_ver, int ind, dpOBJB& alpha_r_z, float& lookup_ind ) {
  // r can be -9999 when find_conv() fails
  if ( r<0 || r>fMaxR ) { lookup_ind = -9999.; return 0; }
  if ( alpha<0 || z_DC_m_z_ver<0 ) {
    cout<<"input error!!!"<<endl;
    cout<<"lookup table contain only e- entries, please convert e+ to e- before fitting with lookup table"<<endl;
    lookup_ind = -9999;
    return 0;
  }
  float lookup_alpha_delta = fLookupAlphaDeltaBase;
  float lookup_r_delta = fLookupRDeltaBase;
  float lookup_z_delta = fLookupZDeltaBase;
  vector<dpOBJA> points;
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
    vector<dpOBJA> points_temp = points;
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


bool dpReco::ProjectPhi( float alpha, float phiDC, float r, float z_DC_m_z_ver, dpOBJB& alpha_r_z, float& phi_r ) {
  //DC coord sys: phi -0.5Pi~1.5Pi
  if ( (r<0) || (r>fMaxR) ) { phi_r=-9999.; return 0; }
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

bool dpReco::ProjectTheta( float alpha, float r, float z_DC_m_z_ver,
				    dpOBJB& alpha_r_z, float& theta_r ) {
  //DC coord sys: phi -0.5Pi~1.5Pi
  if ( (r<0) || (r>fMaxR) ) {
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

bool dpReco::DeltaPhi( float alpha_e, float alpha_p,
				float phi_e, float phi_p,
				float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
				dpOBJB& alpha_r_z, float& dphi ) {
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

bool dpReco::FindR( float alpha_e, float alpha_p,
			    float phi_e, float phi_p,
			    float z_DC_m_z_ver_e, float z_DC_m_z_ver_p,
			    dpOBJB& alpha_r_z, float& r_conv ) {
  // assumption: DeltaPhi(r) should be monotomic increase within [0, fMaxR]
  // hence |DeltaPhi(r)| should be unimodal
  // use godel section search to find the minimum of |DeltaPhi(r)| within range [0, fMaxR]
  double a = 0, b = fMaxR;
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
  // can't do LookupFit for the whole r range [0, fMaxR], dphi always -9999., and that's the minimum value
  return 1;
}

void dpReco::findIntersection(CNTE *pos,
				      CNTE *neg,
				      CNTDE *de,
				      float zvertex) {
  if( pos->GetAlpha()*neg->GetAlpha()>0 )
    return;

  float r_conv, phi_conv_e, phi_conv_p, theta_conv_e, theta_conv_p;
  FindR( pos->GetAlpha(), neg->GetAlpha(),
	 pos->GetPhi(), neg->GetPhi(),
	 pos->GetZed() - zvertex, neg->GetZed() - zvertex,
	 fLUT->alpha_r_z, r_conv );
  ProjectPhi( pos->GetAlpha(), pos->GetPhi(), r_conv, pos->GetZed() - zvertex, fLUT->alpha_r_z, phi_conv_e);
  ProjectPhi( neg->GetAlpha(), neg->GetPhi(), r_conv, neg->GetZed() - zvertex, fLUT->alpha_r_z, phi_conv_p );
  ProjectTheta( pos->GetAlpha(), r_conv, pos->GetZed() - zvertex, fLUT->alpha_r_z, theta_conv_e );
  ProjectTheta( neg->GetAlpha(), r_conv, neg->GetZed() - zvertex, fLUT->alpha_r_z, theta_conv_p );

  de->SetR( r_conv );
  if ( r_conv<-9000. ) {
    de->SetEphi( -9999. );
    de->SetPphi( -9999. );
    de->SetEtheta( -9999. );
    de->SetPtheta( -9999. );
    return;
  }
  de->SetEphi( phi_conv_e );
  de->SetPphi( phi_conv_p );
  de->SetEtheta( theta_conv_e );
  de->SetPtheta( theta_conv_p );
}

TVector3 dpReco::findMomentum(CNTE* trk, float r, float phi_conv, float theta_conv, float zvertex) {
  TVector3 vec;
  float mom, phi, theta;
  LookupFit( TMath::Abs(trk->GetAlpha()), r, TMath::Abs(trk->GetZed() - zvertex), 1, fLUT->alpha_r_z, mom );
  ProjectPhi( trk->GetAlpha(), trk->GetPhi(), r, TMath::Abs(trk->GetZed() - zvertex), fLUT->alpha_r_z, phi );
  ProjectTheta( trk->GetAlpha(), r, trk->GetZed() - zvertex, fLUT->alpha_r_z, theta );
  float pt = mom*sin(theta);
  vec.SetPtThetaPhi( pt, theta, phi );
  return vec;
}
