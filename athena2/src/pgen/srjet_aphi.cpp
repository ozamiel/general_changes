//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file srjet.cpp
//! \brief Sets up a relativistic jet introduced through L-x1 boundary (left edge)
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt(), pi
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"


#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// BCs on L-x3 (lower edge) of grid with jet inflow conditions
void JetInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh);
int RefinementCondition(MeshBlock *pmb);



namespace {
// Make radius of jet and jet variables global so they can be accessed by BC functions
// Real r_amb,
  Real d_amb, p_amb, vx_amb, vy_amb, vz_amb, bx_amb, by_amb, bz_amb;
  Real r_jet, a, d_jet, p_jet, vx_jet, vy_jet, vz_jet, bx_jet, by_jet, bz_jet, b_0, dr_jet;
  Real mang, dang;
  Real gad, gam_add, gm1, x1_0, x2_0, x1min;
  Real atw_jet, atw_amb, hg_amb, hg_jet, rang_jet, rang_amb, phang_jet, phang_amb, d;
  Real SmoothStep(Real x);
  Real A2_intg(Real x1);
  Real A2(Real x1, Real x3);
  Real aintg1(Real x1);
  Real aintg2(Real x1);
    
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
    // ambient medium parameters:
  d_amb  = pin->GetReal("problem", "d");
  p_amb  = pin->GetReal("problem", "p");
  vx_amb = pin->GetReal("problem", "vx");
  vy_amb = pin->GetReal("problem", "vy");
  vz_amb = pin->GetReal("problem", "vz");
  if (MAGNETIC_FIELDS_ENABLED) {// supports uniform ambient MF
    bx_amb = pin->GetReal("problem", "bx");
    by_amb = pin->GetReal("problem", "by");
    bz_amb = pin->GetReal("problem", "bz");
  }
  // inside the jet:
  d_jet  = pin->GetReal("problem", "djet");
  p_jet  = pin->GetReal("problem", "pjet");
  vx_jet = pin->GetReal("problem", "vxjet"); // sets the opening angle of the jet (\tg = vxjet/vyjet)
  vy_jet = pin->GetReal("problem", "vyjet"); // sets the rotation 4-velocity at the jet boundary
  vz_jet = pin->GetReal("problem", "vzjet");
  
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_jet = pin->GetReal("problem", "bxjet");
    by_jet = pin->GetReal("problem", "byjet");
    bz_jet = pin->GetReal("problem", "bzjet");
    b_0  = pin->GetReal("problem", "b0");
  }
  r_jet = pin->GetReal("problem", "rjet");
  dr_jet = pin->GetReal("problem", "drjet");
  x1min = mesh_size.x1min;
  x1_0 = 0.5*(mesh_size.x1max + mesh_size.x1min);
  x2_0 = 0.5*(mesh_size.x2max + mesh_size.x2min);
  
  // openangle = pin->GetReal("problem", "openangle"); // opening angle of the jet, radians
  // angular perturbations of the jet boundary shape
  mang = pin->GetReal("problem", "mang");
  dang = pin->GetReal("problem", "dang");
  
  gad = pin->GetReal("hydro", "gamma"); // adiabatic index
  gam_add = gad/(gad-1.);
  
  // parameter combinations
  Real gamma_amb = sqrt(1.+vx_amb*vx_amb+vy_amb*vy_amb+vz_amb*vz_amb); //ambient Lorentz factor
  atw_amb = gamma_amb*gamma_amb * (d_amb + gad/(gad-1.) * p_amb) ; // Atwood parameter
  hg_amb = (1.+gad/(gad-1.)*p_amb / d_amb) * gamma_amb; // \gamma_inf
  rang_amb = vx_amb / vz_amb ;
  phang_amb = vy_amb / vz_amb ;
  //  Lorentz factor inside the jet
  Real gamma_jet = sqrt(1.+vx_jet*vx_jet+vy_jet*vy_jet+vz_jet*vz_jet);
  atw_jet = gamma_jet*gamma_jet * (d_jet + gad/(gad-1.) * p_jet) ; // Atwood parameter
  hg_jet = (1.+gad/(gad-1.)*p_jet / d_jet) * gamma_jet; // \gamma_inf
  rang_jet = vx_jet / vz_jet ;
  phang_jet = vy_jet / vz_jet ;
  a = r_jet / 2.;
  d = 1./(4.*dr_jet*dr_jet*dr_jet);
  
  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, JetInnerX3);
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Jet problem
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;
  gad = peos->GetGamma() ;
  
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  
  // initialize conserved variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
	phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = d_amb;
	if(std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
	  phydro->w(IM1,k,j,i) = phydro->w1(IM1,k,j,i) = vx_amb;
	  phydro->w(IM2,k,j,i) = phydro->w1(IM2,k,j,i) = vy_amb;
	  phydro->w(IM3,k,j,i) = phydro->w1(IM3,k,j,i) = vz_amb;
	  
	}
	else{
	  phydro->w(IM1,k,j,i) = phydro->w1(IM1,k,j,i) = vx_amb;
	  phydro->w(IM2,k,j,i) = phydro->w1(IM2,k,j,i) = vy_amb; // perpendicular
	  phydro->w(IM3,k,j,i) = phydro->w1(IM3,k,j,i) = vz_amb; // along the jet
	}
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = p_amb;
      }
    }
  }
  
  //for (int k=kl; k<=ku; ++k) {
  //  for (int i=il; i<=iu; ++i) {
      
  //    Real r = pcoord->x1f(i);
  //    Real z = pcoord->x3f(k);
  //    std::cout << r << " " << z << " " << A2(r,z) << "\n";
  //  }
  //}
  
  // AthenaArray<Real> bb;
  
  // initialize interface B
  if (MAGNETIC_FIELDS_ENABLED) {
      
    AthenaArray<Real> area, len, len_p1;
    area.NewAthenaArray(ncells1);
    len.NewAthenaArray(ncells1);
    len_p1.NewAthenaArray(ncells1);
      
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
	pcoord->Face1Area(k,j,il,iu+1,area);
	pcoord->Edge2Length(k  ,j,il,iu+1,len);
	pcoord->Edge2Length(k+1,j,il,iu+1,len_p1);
	for (int i=il; i<=iu+1; ++i) {
	  Real rf = pcoord->x1f(i);
	  Real zf = pcoord->x3f(k);
	  Real zf_p1 = pcoord->x3f(k+1);
	  pfield->b.x1f(k,j,i) = -(len_p1(i)*A2(rf,zf_p1) - len(i)*A2(rf,zf))/area(i);
          
	}
      }
    }
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
	for (int i=il; i<=iu; ++i) {
	  pfield->b.x2f(k,j,i) = by_amb;
	  //  bb(IB2, k,j,i) = by_amb;
	}
      }
    }
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face3Area(k,j,il,iu,area);
        pcoord->Edge2Length(k,j,il,iu+1,len);
        for (int i=il; i<=iu; ++i) {
	  Real rf = pcoord->x1f(i);
	  Real rf_p1 = pcoord->x1f(i+1);
	  Real zf = pcoord->x3f(k);
          pfield->b.x3f(k,j,i) = (len(i+1)*A2(rf_p1,zf) - len(i)*A2(rf,zf))/area(i);
        }
      }
    }
    
    
    
    
    // Calculate cell-centered magnetic field
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, il, iu, jl, ju, kl,
                                       ku);
    
    
    // Initialize conserved values/
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl, ju,
			       kl, ku);
    
    
    // peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void JetInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for jet problem

void JetInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables in inlet ghost zones
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
	Real rad, pert;
	prim(IPR,kl-k,j,i) = p_amb;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          rad = pco->x1v(i);
	  pert = (1.+dang * cos(pco->x2v(j)*mang)) ; // perturbation   
          rad *= pert ; 
        }
        else{
	  rad = std::sqrt(SQR(pco->x1v(i)-x1min) + SQR(pco->x2v(j)-x2_0));
	  //phi = std::atan2(y,x)
	}
	Real step = SmoothStep((rad - r_jet)/dr_jet);
	// Real divfactor = (rad/pert-x1min) / r_jet * openangle ;  // opening * (R/Rjet)
	Real atw = (atw_jet-atw_amb) * step + atw_amb ;
	Real hg = (hg_jet - hg_amb) * step + hg_amb ;
	Real rang = (rang_jet - rang_amb) * step + rang_amb ;
	rang *= (rad - x1min) / r_jet ; // expansion rate -> physical velocity
	Real phang = (phang_jet - phang_amb) * step + phang_amb ;
	phang *= (rad-x1min) / r_jet ; // omega -> physical velocity
	Real R = pco->x1v(i);
	Real smfnc = SmoothStep((R - r_jet)/dr_jet);
	
        
        Real p = prim(IPR,kl-k,j,i);
	
	Real gamma_amb = sqrt(1.+vx_amb*vx_amb+vy_amb*vy_amb+vz_amb*vz_amb);
	Real gamma_jet = sqrt(1.+vx_jet*vx_jet+vy_jet*vy_jet+vz_jet*vz_jet);
	Real Bphi_const = (b_0*a*r_jet/(a*a+r_jet*r_jet)) * smfnc;
	
	Real smfnc_c = SmoothStep((x1min - r_jet)/dr_jet);
	Real b_phi_cen = (b_0*a*x1min/(a*a+x1min*x1min)) * smfnc_c;
	Real p_cen = p_amb; //Setting this because prim(IPR,kl-k,0,0) would bring not the pressure at the center but at the beginning of the domain
	Real bern_jet = (1. + (gam_add*p_cen/d_jet))*gamma_jet + (1/(gamma_jet*d_jet))*(b_phi_cen*b_phi_cen); //Setting the Bernoulli parameter and smoothing it
	Real bern_amb = (1. + (gam_add*p_amb/d_amb))*gamma_amb;
	Real bern_sm = (bern_jet - bern_amb) * smfnc + bern_amb;
	
	Real atwd_jet = gamma_jet*gamma_jet*(d_jet + gam_add*p_cen); //Setting the Atwood parameter and smoothing it
	Real atwd_amb = gamma_amb*gamma_amb*(d_amb + gam_add*p_amb);
	Real atwd_sm = (atwd_jet - atwd_amb) * smfnc + atwd_amb;
	
	Real Psi = (atwd_sm + Bphi_const*Bphi_const)/bern_sm; //Setting a parameter Psi to calculate gamma and rho from it
	
	Real gamma = (Psi/(2.*gam_add*p)) * (sqrt(1. + (4.*gam_add*p*atwd_sm)/(Psi*Psi)) - 1.);
	prim(IVZ,kl-k,j,i) = sqrt((gamma*gamma-1.)/(1. + rang * rang + phang * phang)); // gamma^2 - 1 = uz^2 + uy^2 + ux^2 
	prim(IVX,kl-k,j,i) = prim(IVZ,kl-k,j,i) * rang; // rang = (rang_jet - rang_amb) * step + rang_amb , rang_i = vx_i/vz_i
	prim(IVY,kl-k,j,i) = prim(IVZ,kl-k,j,i) * phang; // phang = (phang_jet - phang_amb) * step + phang_amb , phang_i = vy_i/vz_i
	prim(IDN,kl-k,j,i) = Psi/gamma;
	
	
      }
    }
  }
  
  if(MAGNETIC_FIELDS_ENABLED) {
    
    Real ncells1 = pmb->ncells1;
    AthenaArray<Real> area, len, len_p1;
    area.NewAthenaArray(ncells1);
    len.NewAthenaArray(ncells1);
    len_p1.NewAthenaArray(ncells1);
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
	pco->Face1Area(k,j,il,iu+1,area);
	pco->Edge2Length(k  ,j,il,iu+1,len);
	pco->Edge2Length(k+1,j,il,iu+1,len_p1);
	for (int i=il; i<=iu+1; ++i) {
	  Real rf = pco->x1f(i);
	  Real zf = pco->x3f(k);
	  Real zf_p1 = pco->x3f(k+1);
	  b.x1f(kl-k,j,i) = -(len_p1(i)*A2(rf,zf_p1) - len(i)*A2(rf,zf))/area(i);
          
	}
      }
    } 
    
    
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
	for (int i=il; i<=iu; ++i) {
	  Real R = pco->x1f(i);
	  Real smfnc = SmoothStep((R - r_jet)/dr_jet);
	  b.x2f(kl-k,j,i) = (b_0*a*R/(a*a+R*R)) * smfnc; //Setting Bphi = R/a * Bz - force-free solution (total pressure-free in our case).
	}
      }
    }
    
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        pco->Face3Area(k,j,il,iu,area);
        pco->Edge2Length(k,j,il,iu+1,len);
        for (int i=il; i<=iu; ++i) {
	  Real rf = pco->x1f(i);
	  Real rf_p1 = pco->x1f(i+1);
	  Real zf = pco->x3f(k);
          b.x3f(kl-k,j,i) = (len(i+1)*A2(rf_p1,zf) - len(i)*A2(rf,zf))/area(i);
        }
      }
    }
    
  }
}

namespace {
  Real SmoothStep(Real x)
  {
    // step function approximation
    
      Real modx = std::max(std::min(x,1.),-1.);
      
    return 1./2. - modx*(3.-modx*modx)/4.;
  }
  
  Real A2_intg(Real x1)
  {
    // an analytic function of the vector potential Aphi divided into 3 section: rmin < r < r_jet - dr_jet, r_jet - dr_jet < r < r_jet + dr_jet, r > r_jet + dr_jet
    Real aphir;
    if(x1<r_jet - dr_jet){
      
      aphir = aintg1(x1) - aintg1(x1min);
    } else if(x1>=r_jet - dr_jet && x1<r_jet + dr_jet){
      aphir = aintg2(x1) - aintg2(r_jet - dr_jet) + aintg1(r_jet - dr_jet) - aintg1(x1min);
    } else {
      aphir = aintg2(r_jet + dr_jet) - aintg2(r_jet - dr_jet) + aintg1(r_jet - dr_jet) - aintg1(x1min);
    }
      
    
    return aphir;
  }
  
  Real A2(Real x1, Real x3) //calcualtes the value of Aphi based on location
  {
    Real r0 = (4.*r_jet*x1 + x1min*x3)/(x3 + 4.*r_jet);
    return A2_intg(r0)/x1;
  }
  
  Real aintg1(Real x1) //aphi for the first section
  {
    return b_0*((a*a)/2.)*log(a*a + x1*x1);
  }
  
  Real aintg2(Real x1) //aphi for the second section
  {
    return b_0*((d*a*a)/6.)*(x1*(-6.*a*a - 18.*dr_jet*dr_jet + 18.*r_jet*r_jet - 9.*r_jet*x1 + 2.*x1*x1) + 6.*a*(a*a + 3.*dr_jet*dr_jet - 3.*r_jet*r_jet)*atan(x1/a) + (9.*r_jet*a*a + 6*dr_jet*dr_jet*dr_jet + 9.*r_jet*dr_jet*dr_jet - 3.*r_jet*r_jet*r_jet)*log(a*a + x1*x1));
  }
  
} // namespace

//Defining a refinement condition - takes the differences of rho and refines the cells in which the difference is greater than 20%, derefines if less than 5%. 
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> &bcc = pmb->pfield->bcc;
  Real maxsig=0.0;
  for(int k=pmb->ks; k<=pmb->ke;k++){
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
	Real Bsqc = bcc(IB1,k,j,i)*bcc(IB1,k,j,i) + bcc(IB2,k,j,i)*bcc(IB2,k,j,i) + bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
	Real sigma_M = Bsqc/w(IDN,k,j,i);
	//Real dif_r = (std::abs(w(IDN,k,j,i+1)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x1v(i+1)-pco->x1v(i));
	//Real dif_z = (std::abs(w(IDN,k+1,j,i)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x3v(k+1)-pco->x3v(k));
	//Real dif_phi = (std::abs(w(IDN,k,j+1,i)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x1v(i)*(pco->x2v(j+1)-pco->x2v(j)));
	if(maxsig < sigma_M) maxsig = sigma_M;
	//if(maxdif < dif_z) maxdif = dif_z;
	//if(maxdif < dif_phi) maxdif = dif_phi;
	
      }
    }
  }
  if(maxsig > 0.01) return 1;
  //if(maxsig < 0.001) return -1;
  return 0;
}
