# dp

variable start with "f " means data member(all variable down have a f on them)

CNTDE central tracker dielectron object
 > Ephi,Pphi, Etheta, Pttheta, R, Eid, Pid

CNTE central tracker electron
 > DCarm,DCside,Px,Py,Pz,Phi,Phi0,Theta0,Zed,Alpha,Charge,EMCid,Ecore,
 Disp,N0,Chi2,Npe0,Dep,Prob,EMCdz,EMCdphi,EMCx,EMcy,EMCz,CenterPhi,CenterZ,
 PPC1x,PPC1y,PPC1z,Esect,Ysect,Zsect

EMCC electromagnetic calorimeter
> Arm,ID,X,Y,Z,Ecore,ProbPhoton,EMCdz,EMCdphi

PhotonEvent
> BBCcharge,ZVertex,Centrality,Psi_BBC,Psi_BBN,Psi_BBCS,Svxz,Run,
 vector<CNTE> P, vector<CNTE> N,vector<CNTDE> DE, vector<EMCC> EMCC
