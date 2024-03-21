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
#include <cmath>      // sqrt()
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
  Real r_jet, a, d_jet, p_jet, vx_jet, vy_jet, vz_jet, bx_jet, by_jet, bz_jet, b_0;
  Real dr_jet;
  Real mang, dang;
  Real gad, gm1, x1_0, x2_0, x1min;
  Real atw_jet, atw_amb, hg_amb, hg_jet, rang_jet, rang_amb, phang_jet, phang_amb;
  Real SmoothStep(Real x);
  Real integ_A3(Real x1_m, Real x1_p);
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
  b_0  = pin->GetReal("problem", "b0");
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
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
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

  // AthenaArray<Real> bb;
    
  // initialize interface B
    if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie+1; ++i) {
	    pfield->b.x1f(k,j,i) = bx_amb;
	    //bb(IB1, k,j,i) = bx_amb;
	  }
	}
      }
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je+1; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    pfield->b.x2f(k,j,i) = by_amb;
	    //  bb(IB2, k,j,i) = by_amb;
	  }
	}
      }
      for (int k=ks; k<=ke+1; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    pfield->b.x3f(k,j,i) = bz_amb;
	    //  bb(IB3, k,j,i) = bz_amb;
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
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          rad = pco->x1v(i);
	  pert = (1.+dang * cos(pco->x2v(j)*mang)) ; // perturbation   
          rad *= pert ; 
        }
        else{
              rad = std::sqrt(SQR(pco->x1v(i)-x1min) + SQR(pco->x2v(j)-x2_0));
	      //phi = std::atan2(y,x)
	}
	Real step = SmoothStep(-(rad - r_jet)/dr_jet);
	// Real divfactor = (rad/pert-x1min) / r_jet * openangle ;  // opening * (R/Rjet)
	Real atw = (atw_jet-atw_amb) * step + atw_amb ;
	Real hg = (hg_jet - hg_amb) * step + hg_amb ;
	Real rang = (rang_jet - rang_amb) * step + rang_amb ;
	rang *= (rad - x1min) / r_jet ; // expansion rate -> physical velocity
	Real phang = (phang_jet - phang_amb) * step + phang_amb ;
	phang *= (rad-x1min) / r_jet ; // omega -> physical velocity
	//Real press = (p_jet - p_amb) * step + p_amb ;
	//Real gamma = atw / (8.*press * hg) * (sqrt(1.+16.*press*hg*hg/atw)-1.);
	//prim(IDN,kl-k,j,i) = atw/gamma/gamma-4.*press ; // (d_jet-d_amb) * step + d_amb;
	//prim(IVZ,kl-k,j,i) =  sqrt((gamma*gamma-1.)/(1. + rang * rang + phang * phang)) ; // (vz_jet-vz_amb) * step + vz_amb;
	//prim(IVX,kl-k,j,i) = prim(IVZ,kl-k,j,i) * rang ; // (vz_jet*divfactor-vx_amb) * step + vx_amb; // radial velocity
	//prim(IVY,kl-k,j,i) = prim(IVZ,kl-k,j,i) * phang ;  // (vy_jet*rad-vy_amb) * step + vy_amb;
	//prim(IPR,kl-k,j,i) = (p_jet-p_amb) * step + p_amb;
	  
		    // prim(IVX,kl-k,j,i) = (vz_jet*divfactor-vx_amb) * step + vx_amb; // radial velocity
		    // prim(IVY,kl-k,j,i) = (vy_jet*rad-vy_amb) * step + vy_amb;
		    // prim(IVZ,kl-k,j,i) = (vz_jet-vz_amb) * step + vz_amb;
		    // prim(IPR,kl-k,j,i) = (p_jet-p_amb) * step + p_amb;
    //Real R = pco->x1v(i);
    Real smfnc = SmoothStep((r_jet - R)/dr_jet);

    Real gamma_amb = sqrt(1.+vx_amb*vx_amb+vy_amb*vy_amb+vz_amb*vz_amb);
    Real gamma_jet = sqrt(1.+vx_jet*vx_jet+vy_jet*vy_jet+vz_jet*vz_jet);

    Real smfnc_c = SmoothStep((r_jet - x1min)/dr_jet);
    Real b_cen = (b_0*x1min/a) * smfnc_c;
    Real p_cen = p_amb - b_cen*b_cen/2; //Setting this because prim(IPR,kl-k,0,0) would bring not the pressure at the center but at the beginning of the domain
    //prim(IPR,kl-k,j,i) = p_amb - b.x2f(kl-k,j,i)*b.x2f(kl-k,j,i)/(2.)
    prim(IPR,kl-k,j,i) = (p_cen-p_amb) * step + p_amb;
    Real press = prim(IPR,kl-k,j,i);
    Real bern_jet = (1. + ((gad/(gad-1.))*p_cen/d_jet))*gamma_jet - b_cen*b_cen/(d_jet*vz_jet); //Setting the Bernoulli parameter and smoothning it
    Real bern_amb = (1. + ((gad/(gad-1.))*p_amb/d_amb))*gamma_amb;
    Real bern_sm = (bern_jet - bern_amb) * smfnc + bern_amb;
          
    Real atwd_jet = gamma_jet*gamma_jet*(d_jet + (gad/(gad-1.))*p_cen); //Setting the Atwood parameter and smoothning it
    Real atwd_amb = gamma_amb*gamma_amb*(d_amb + (gad/(gad-1.))*p_amb);
    Real atwd_sm = (atwd_jet - atwd_amb) * smfnc + atwd_amb;

	//Real Psi = (atwd_sm - b.x2f(kl-k,j,i)*b.x2f(kl-k,j,i))/bern_sm; //Setting a parameter Psi to calculate gamma and rho from it
    Real Psi = (atwd_sm)/bern_sm;

    Real gamma = (Psi/(2.*(gad/(gad-1.))*press)) * (sqrt(1. + (4.*(gad/(gad-1.))*press*atwd_sm)/(Psi*Psi)) - 1.);
    prim(IVZ,kl-k,j,i) = sqrt((gamma*gamma-1.)/(1. + rang * rang + phang * phang)); // gamma^2 - 1 = uz^2 + uy^2 + ux^2 
    prim(IVX,kl-k,j,i) = prim(IVZ,kl-k,j,i) * rang; // rang = (rang_jet - rang_amb) * step + rang_amb , rang_i = vx_i/vz_i
    prim(IVY,kl-k,j,i) = prim(IVZ,kl-k,j,i) * phang; // phang = (phang_jet - phang_amb) * step + phang_amb , phang_i = vy_i/vz_i
    prim(IDN,kl-k,j,i) = Psi/gamma;
    
      }
    }
  }
}

namespace {
  Real SmoothStep(Real x)
  {
    // step function approximation
    return (tanh(x)+1.)/2. ; // x/std::sqrt(x*x+1.);
  }
  
  Real integ_A3(Real x1_m, Real x1_p)
  {
    // calculating the difference of the vector potential between x1_m and x1_p.
    Real tot_steps = 100;
    Real dr = (x1_p - x1_m)/tot_steps;
    Real integ = 0.;
    for(int step_num=1 ; step_num<=tot_steps ; ++step_num) {
      Real r = x1_m + step_num*dr;
      if(x1_p<a) {
	integ += (b_0/(2.*a)) * (tanh((r_jet - r)/dr_jet) + 1.) * r * dr;
      } else {
	integ += (b_0/(2.*r)) * (tanh((r_jet - r)/dr_jet) + 1.) * a * dr;
      }
    }
    return integ;
  }
} // namespace

//Defining a refinement condition - takes the differences of rho and refines the cells in which the difference is greater than 20%, derefines if less than 5%. 
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxdif=0.0;
  for(int k=pmb->ks; k<=pmb->ke;k++){
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
	Real dif_r = (std::abs(w(IDN,k,j,i+1)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x1v(i+1)-pco->x1v(i));
	Real dif_z = (std::abs(w(IDN,k+1,j,i)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x3v(k+1)-pco->x3v(k));
	Real dif_phi = (std::abs(w(IDN,k,j+1,i)-w(IDN,k,j,i)))/w(IDN,k,j,i);   //(pco->x1v(i)*(pco->x2v(j+1)-pco->x2v(j)));
	//double list[] = {maxdif, dif_r, dif_z, dif_phi};
	//if maxdif = std::max_element(list,list+4);
	if(maxdif < dif_r) maxdif = dif_r;
	if(maxdif < dif_z) maxdif = dif_z;
	if(maxdif < dif_phi) maxdif = dif_phi;

      }
    }
  }
  if(maxdif > 0.2) return 1;
  if(maxdif < 0.05) return -1;
  return 0;
}
