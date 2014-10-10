//
// Reads CAMB matter power spectrum camb_matterpower.dat, and
// provides power spectrum to 2LPT calculation in lpt.c
//
// Based on N-GenIC power.c by Volker Springel
//   http://www.mpa-garching.mpg.de/gadget/right.html#ICcode
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "common.h"

#define WORKSIZE 100000

static double Omega, OmegaLambda;

#define N_A_GROWTH 1000
static struct growth_table
{
  double a,gf;
} *GrowthTable;

static double omega_a(const double a)
{
  return Omega/(Omega+(1-Omega-OmegaLambda)*a+
		 OmegaLambda*a*a*a);
}

double Qfactor(const double a)
{
  double acube=a*a*a;
  return sqrt(Omega*acube+OmegaLambda*acube*acube+
	      (1-Omega-OmegaLambda)*acube*a);
}

double GrowthFactor(const double a)
{
  if(a<0||a>1)
    msg_abort(3000,"Wrong scale factor %lf\n",a);
  if(a==1) return 1;
  else if(a==0) return 0;
  else {
    int ibin=(int)(a*N_A_GROWTH);
    double u=(a-GrowthTable[ibin].a)*N_A_GROWTH;
    double gfac=(1-u)*GrowthTable[ibin].gf+u*GrowthTable[ibin+1].gf;

    return gfac;
  }
}

double GrowthFactor2(const double a)
{
  double gf1=GrowthFactor(a);
  double om=omega_a(a);

  return -0.42857142857*gf1*gf1*pow(om,-0.00699300699); //-3/7*D_1^2*Omega_M^{-1/143}
}

double Vgrowth(const double a) //Returns dD1/dTau
{
  double om=omega_a(a);
  double q=Qfactor(a);
  double f1=pow(om,0.555);
  double d1=GrowthFactor(a);
  
  return q*d1*f1/a;
}

double Vgrowth2(const double a) //Returns dD2/dTau
{
  double om=omega_a(a);
  double q=Qfactor(a);
  double f2=2*pow(om,0.545454);
  double d2=GrowthFactor2(a);
  
  return q*d2*f2/a;
}

static double growth_int(double a,void *param)
{
  return pow(a/(Omega+(1-Omega-OmegaLambda)*a+OmegaLambda*a*a*a),1.5);
}

static double growth(double a)
{
  double hubble_a;

  hubble_a=sqrt(Omega/(a*a*a)+(1-Omega-OmegaLambda)/(a*a)+OmegaLambda);

  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace=gsl_integration_workspace_alloc(WORKSIZE);

  F.function=&growth_int;

  gsl_integration_qag(&F,0,a,0,1.0e-8,WORKSIZE,GSL_INTEG_GAUSS41, 
		      workspace,&result,&abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a*result;
}

static void get_growth_factor(void)
{
  int ii;
  double g0=growth(1.0);
  GrowthTable=malloc(N_A_GROWTH*sizeof(struct growth_table));

  for(ii=0;ii<N_A_GROWTH;ii++) {
    double a=(double)ii/N_A_GROWTH;
    double gf=growth(a)/g0;

    GrowthTable[ii].a=a;
    GrowthTable[ii].gf=gf;
  }
}

void cosmo_init(const double omega_m,const double omega_lambda)
{
  Omega=omega_m;
  OmegaLambda=omega_lambda;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  if(myrank==0)
    get_growth_factor();

  if(myrank!=0)
    GrowthTable=malloc(N_A_GROWTH*sizeof(struct growth_table));

  MPI_Bcast(GrowthTable,N_A_GROWTH*sizeof(struct growth_table),MPI_BYTE,0,
	    MPI_COMM_WORLD);
}
