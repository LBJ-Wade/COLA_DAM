#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>

#include "common.h"

Particles *allocate_particles(const unsigned long long np_total,
			      const double np_alloc_factor)
{
  Particles *particles=malloc(sizeof(Particles));
  const int nnode=comm_nnode();
  const int np_alloc=(int)((double)(np_alloc_factor*np_total)/nnode);

  particles->p=malloc(sizeof(Particle)*np_alloc);
  if(particles->p==0)
    msg_abort(0010,"Error: Failed to allocate memory for particles\n");

  particles->force=calloc(sizeof(float),3*np_alloc);
  if(particles->force==0)
    msg_abort(0010,"Error: Failed to allocate memory for particle forces\n");

  msg_printf(info,"%d Mbytes allocated for %d particles (alloc_factor= %.2lf)\n",
	     (sizeof(Particle)+3*sizeof(float))*np_alloc/(1024*1024),
	     np_alloc,np_alloc_factor);

  particles->np_allocated=np_alloc;

  particles->np_total=np_total;
  particles->np_average=(float)((double)np_total/nnode);

  return particles;
}

Snapshot *allocate_snapshot(const unsigned long long np_total,const int np_alloc,
			    void* const mem,const size_t mem_size)
{
  Snapshot* snapshot=malloc(sizeof(Snapshot));

  snapshot->np_allocated=np_alloc;
  snapshot->np_total=np_total;
  snapshot->p=mem; assert(mem_size>=sizeof(ParticleMinimum)*np_alloc);
  snapshot->a=0.0f; //snapshot->a_v= 0.0f; snapshot->a_x= 0.0f;

  return snapshot;
}

void allocate_shared_memory(const int n_grid,Memory* const mem)
{
  // Allocate shared memory

  // mem1
  // Memory for PM (nc_factor^3 * np_local each)
  const int Ngrid=n_grid;
  ptrdiff_t local_nx,local_x_start;
  ptrdiff_t local_ny,local_y_start;
  ptrdiff_t size_pm_one=
    fftwf_mpi_local_size_3d_transposed(Ngrid,Ngrid,Ngrid,MPI_COMM_WORLD,
				       &local_nx,&local_x_start,
				       &local_ny,&local_y_start);
#ifdef _DAM_SAVEMEM
  size_pm_one/=Ngrid;
  size_pm_one*=(Ngrid/2+1);
#endif //_DAM_SAVEMEM
  ptrdiff_t ncomplex_pm=size_pm_one;

  msg_printf(verbose,"%d Mbytes requested for one PM grid\n",
	     (int)(ncomplex_pm*sizeof(fftwf_complex)/(1024*1024)));

  msg_printf(verbose,"PM size %d %d %d\n",size_pm_one,
	     local_nx*Ngrid*Ngrid,local_nx);

  ptrdiff_t ncomplex1=ncomplex_pm;

  mem->mem1=fftwf_alloc_complex(ncomplex1);
  mem->size1=sizeof(fftwf_complex)*ncomplex1;

  if(mem->mem1==0)
    msg_abort(0050,"Error: Unable to allocate %d Mbytes for mem1\n",
	      (int)(mem->size1/(1024*1024)));

  // mem2
  // PM density_k mesh and snapshot
  size_t ncomplex2=(Ngrid/2+1)*Ngrid*local_ny; //ncomplex_pm;
  size_t size2=sizeof(fftwf_complex)*(Ngrid/2+1)*Ngrid*local_ny;

  msg_printf(verbose,"%d Mbytes requested for delta_k mesh (mem2). ny=%d\n",
	     (int)(size2/(1024*1024)),local_ny);

  mem->mem2=fftwf_alloc_complex(ncomplex2);
  mem->size2=sizeof(fftwf_complex)*ncomplex2;

  if(mem->mem2==0)
    msg_abort(0060,"Error: Unable to allocate %d + %d Mbytes for mem1&2.\n",
	      (int)(mem->size1/(1024*1024)),(int)(mem->size2/(1024*1024)));

  msg_printf(info,"%d Mbytes allocated for mem1.\n", 
	     (int)(mem->size1/(1024*1024)));
  msg_printf(info,"%d Mbytes allocated for mem2.\n", 
	     (int)(mem->size2/(1024*1024)));
}
