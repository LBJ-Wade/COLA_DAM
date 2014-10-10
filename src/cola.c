//
// COLA time integration using given force and 2LPT displacement
//
// This code is a modification to the original serial COLA code
// by Svetlin Tassev. See below.
//

/*
    Copyright (c) 2011-2013       Svetlin Tassev
                           Harvard University, Princeton University
 
    This file is part of COLAcode.

    COLAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    COLAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with COLAcode.  If not, see <http://www.gnu.org/licenses/>.
*/



/*

    This is COLAcode: a serial particle mesh-based N-body code 
     illustrating the COLA (COmoving Lagrangian Acceleration) method 
     described in S. Tassev, M. Zaldarriaga, D. Eisenstein (2012).
     Check that paper (refered to as TZE below) for the details. 
     Before using the code make sure you read the README file as well as
     the Warnings section below.
    
    This version: Dec 18, 2012


*/

#include <math.h>
#include <assert.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include "common.h"

static const float nLPT= -2.5f;
static const int fullT= 1; // velocity growth model

//
static double Sphi(double ai, double af, double aRef);
static double Sq(double ai, double af, double aRef);

// Leap frog time integration
// ** Total momentum adjustment dropped

void cola_kick(Particles* const particles,const float avel1)
{
  timer_start(evolve);  
  const float AI=particles->a_v;  // t - 0.5*dt
  const float A= particles->a_x;  // t
  const float AF=avel1;           // t + 0.5*dt

  msg_printf(normal,"Kick %g -> %g\n",AI,avel1);

  const float Om=Param.omega_m;
  const float dda=Sphi(AI,AF,A);

  const float q1=GrowthFactor(A);
  const float q2=(GrowthFactor2(A)-q1*q1);

  msg_printf(normal,"growth factor %g\n",q1);

  Particle* const P=particles->p;
  const int np=particles->np_local;
  float3* const f=particles->force;
  
  // Kick using acceleration at a= A
  // Assume forces at a=A is in particles->force

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0; i<np; i++) {
    float ax=-1.5*Om*dda*(f[i][0]+P[i].dx1[0]*q1+P[i].dx2[0]*q2);
    float ay=-1.5*Om*dda*(f[i][1]+P[i].dx1[1]*q1+P[i].dx2[1]*q2);
    float az=-1.5*Om*dda*(f[i][2]+P[i].dx1[2]*q1+P[i].dx2[2]*q2);

    P[i].v[0]+=ax;
    P[i].v[1]+=ay;
    P[i].v[2]+=az;
  }

  //velocity is now at a= avel1
  particles->a_v=avel1;
  timer_stop(evolve);  
}

void cola_drift(Particles* const particles,const float apos1)
{
  timer_start(evolve);
  const float A= particles->a_x; // t
  const float AC=particles->a_v; // t + 0.5*dt
  const float AF=apos1;          // t + dt

  Particle* const P=particles->p;
  const int np=particles->np_local;
  
  const float dyyy=Sq(A, AF, AC);

  const float da1=GrowthFactor(AF)-GrowthFactor(A);    // change in D_{1lpt}
  const float da2=GrowthFactor2(AF)-GrowthFactor2(A);  // change in D_{2lpt}

  msg_printf(normal,"Drift %g -> %g\n",A,AF);
    
  // Drift
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0;i<np;i++) {
    P[i].x[0]+=P[i].v[0]*dyyy+P[i].dx1[0]*da1+P[i].dx2[0]*da2;
    P[i].x[1]+=P[i].v[1]*dyyy+P[i].dx1[1]*da1+P[i].dx2[1]*da2;
    P[i].x[2]+=P[i].v[2]*dyyy+P[i].dx1[2]*da1+P[i].dx2[2]*da2;
  }
    
  particles->a_x=AF;
  timer_stop(evolve);
}

//
// Functions for our modified time-stepping (used when StdDA=0):
//
static double gpQ(double a) { 
  return pow(a,nLPT);
}

static double fun (double a,void *params) {
  double f;
  if (fullT==1) f=gpQ(a)/Qfactor(a); 
  else f=1.0/Qfactor(a);
  
  return f;
}

static double Sq(double ai,double af,double aRef) {
  gsl_integration_workspace *w 
    =gsl_integration_workspace_alloc(5000);
  
  double result,error;
  double alpha=0;
  
  gsl_function F;
  F.function=&fun;
  F.params=&alpha;
  
  gsl_integration_qag(&F,ai,af,0,1e-5,5000,6,
		      w,&result,&error); 
  
  gsl_integration_workspace_free(w);
     
  if (fullT==1)
    return result/gpQ(aRef);
  return result;
}
     
static double DERgpQ(double a) { // This must return d(gpQ)/da
  return nLPT*pow(a,nLPT-1);
}
     
static double Sphi(double ai,double af,double aRef) {
  double result;
  result=(gpQ(af)-gpQ(ai))*aRef/Qfactor(aRef)/DERgpQ(aRef);
  
  return result;
}

// Interpolate position and velocity for snapshot at a=aout
void cola_set_snapshot(const double aout,Particles const * const particles,
			      Snapshot* const snapshot)
{
  timer_start(interp);
  const int np=particles->np_local;
  Particle const * const p=particles->p;
  float3* const f=particles->force;

  ParticleMinimum* const po=snapshot->p;
  const float Om=Param.omega_m; assert(Om>=0.0f);

  msg_printf(verbose,"Setting up snapshot at a= %4.2f (z=%4.2f) <- %4.2f %4.2f.\n",
	     aout,1.0f/aout-1,particles->a_x,particles->a_v);

  const float vfac=100.0f/aout;   // km/s; H0= 100 km/s/(h^-1 Mpc)

  const float AI=particles->a_v;
  const float A= particles->a_x;
  const float AF=aout;

  const float dda=Sphi(AI,AF,A);

  msg_printf(normal,"Growth factor of snapshot %f (a=%.3f)\n",GrowthFactor(AF),AF);

  const float q1=GrowthFactor(A);          //DAM: might be better to use 0.5*(av+aout) instead of A
  const float q2=(GrowthFactor2(A)-q1*q1); //DAM: but let's leave it for the moment

  const float Dv=Vgrowth(aout); // dD_{za}/dTau
  const float Dv2=Vgrowth2(aout); // dD_{2lpt}/dTau

  const float AC=particles->a_v;
  const float dyyy=Sq(A, AF, AC); //DAM: I'm not sure about the central value in this
                                  //DAM: and dda, but let's leave like this for the moment

  msg_printf(debug,"velocity factor %e %e\n",vfac*Dv,vfac*Dv2);
  msg_printf(debug, "RSD factor %e\n", aout/Qfactor(aout)/vfac);

  const float da1=GrowthFactor(AF)-GrowthFactor(A);    // change in D_{1lpt}
  const float da2=GrowthFactor2(AF)-GrowthFactor2(A);  // change in D_{2lpt}

#ifdef _OPENMP
  #pragma omp parallel for default(shared)  
#endif
  for(int i=0; i<np; i++) {
    // Kick + adding back 2LPT velocity + convert to km/s
    float ax=-1.5*Om*dda*(f[i][0]+p[i].dx1[0]*q1+p[i].dx2[0]*q2);
    float ay=-1.5*Om*dda*(f[i][1]+p[i].dx1[1]*q1+p[i].dx2[1]*q2);
    float az=-1.5*Om*dda*(f[i][2]+p[i].dx1[2]*q1+p[i].dx2[2]*q2);

    po[i].v[0]=vfac*(p[i].v[0]+ax+p[i].dx1[0]*Dv+p[i].dx2[0]*Dv2);
    po[i].v[1]=vfac*(p[i].v[1]+ay+p[i].dx1[1]*Dv+p[i].dx2[1]*Dv2);
    po[i].v[2]=vfac*(p[i].v[2]+az+p[i].dx1[2]*Dv+p[i].dx2[2]*Dv2);

    // Drift
    po[i].x[0]=p[i].x[0]+p[i].v[0]*dyyy+p[i].dx1[0]*da1+p[i].dx2[0]*da2;
    po[i].x[1]=p[i].x[1]+p[i].v[1]*dyyy+p[i].dx1[1]*da1+p[i].dx2[1]*da2;
    po[i].x[2]=p[i].x[2]+p[i].v[2]*dyyy+p[i].dx1[2]*da1+p[i].dx2[2]*da2;

    po[i].id=p[i].id;
  }

  snapshot->np_local=np;
  snapshot->a=aout;
  timer_stop(interp);
}
