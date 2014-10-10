#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"

//#define LOGTIMESTEP 1

int mpi_init(int* p_argc,char ***p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout,const int iout, 
		   Particles const * const particles, 
		   Snapshot * const snapshot,
		   Memory mem);

int main(int argc,char* argv[])
{
  const int multi_thread=mpi_init(&argc,&argv);
  msg_init();
  timer_set_category(Init);

  //
  // Initialization / Memory allocation
  //						      
  if(argc<2)
    msg_abort(1,"Error: Parameter file not specified. cola_halo param.init\n");
  read_parameters(argv[1]);

  msg_set_loglevel(Param.loglevel);

  fft_init(multi_thread);
  comm_init(Param.n_grid,Param.boxsize);

  cosmo_init(Param.omega_m,1-Param.omega_m);

  const double a_init=Param.a_init;
  const double a_final=Param.aout[Param.n_aout-1];
  const double da=(a_final-a_init)/(Param.ntimestep+0);

  Memory mem; 
  allocate_shared_memory(Param.n_grid,&mem);


  Particles* particles=allocate_particles(Param.n_part,Param.np_alloc_factor);
  Snapshot* snapshot= allocate_snapshot(Param.n_part,particles->np_allocated,
					mem.mem2,mem.size2);
  
  pm_init(Param.n_grid,Param.n_part,Param.boxsize,
	  mem.mem1,mem.size1,mem.mem2,mem.size2);

  for(int i=0;i<Param.n_aout;i++) {
    msg_printf(verbose,"zout[%d]= %lf, aout= %f\n", 
	       i,1.0/Param.aout[i]-1,Param.aout[i]);
    if(Param.aout[i]<a_init) {
      msg_abort(1234,"Output required before initial condition %lE > %lE!\n",
		a_init,Param.aout[i]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  // Sets initial grid and 2LPT displacement
  timer_set_category(Reading);
  read_initial_condition(particles);
  // velocities must be zero
  // set av and ax to a_init
  //
  //  lpt_set_initial_condition(a_init,Param.boxsize,particles); //DAM Read here

  timer_set_category(COLA);
  msg_printf(info,"timestep linear in a\n");

  //
  // Time evolution loop
  //
  int iout=0;
  int n_aout=Param.n_aout;
  if(n_aout>0 && Param.ntimestep>1 && a_final>a_init) {
    msg_printf(normal,"Time integration a= %g -> %g, %d steps\n", 
	       a_init,a_final,Param.ntimestep);
    for (int istep=0;istep<Param.ntimestep;istep++) { 
      double avel0=particles->a_v;
      double apos0=particles->a_x;
      double avel1=a_init+(istep+0.5)*da;
      double apos1=a_init+(istep+1.0)*da;
      
      msg_printf(normal,"Timestep %d/%d\n",istep,Param.ntimestep);
      timer_start(comm);
      // move particles to other nodes
      move_particles2(particles,Param.boxsize,mem.mem1,mem.size1);
      timer_stop(comm);
      pm_calculate_forces(particles);

      if(avel1<=avel0 || apos1<=apos0 || apos1<=avel1)
	msg_abort(2,"Something went wrong...\n");

      while(iout<Param.n_aout && avel0<=Param.aout[iout] && Param.aout[iout]<=apos0) {
	snapshot_time(Param.aout[iout],iout,particles,snapshot,mem);
	iout++;
      }
      if(iout>=Param.n_aout) break;
      if(avel1>a_final) break;

      // Leap-frog "kick" -- velocities updated
      cola_kick(particles,avel1);

      while(iout<Param.n_aout && apos0<Param.aout[iout] && Param.aout[iout]<=avel1) {
	snapshot_time(Param.aout[iout],iout,particles,snapshot,mem);
	iout++;
      }
      if(iout>=Param.n_aout) break;

      // Leap-frog "drift" -- positions updated
      cola_drift(particles,apos1);
    }
  }

  timer_print();

  MPI_Finalize();
  return 0;
}

int mpi_init(int* p_argc,char*** p_argv)
{
  // MPI+OpenMP paralellization: MPI_THREAD_FUNNELED
  // supported by mpich2 1.4.1, but now by openmpi 1.2.8

#ifdef _OPENMP
  int thread_level, hybrid_parallel;
  MPI_Init_thread(p_argc,p_argv,MPI_THREAD_FUNNELED,&thread_level);
  hybrid_parallel=(thread_level>=MPI_THREAD_FUNNELED);

  int myrank; MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  if(myrank==0) {
    if(hybrid_parallel)
      printf("MPI + multi thread supported (MPI_THREAD_FUNNELED).\n");
    else
      printf("Warning: MPI + multi thread not supported. 1 thread per node.\n");
  }
	
  return hybrid_parallel;
#else
  MPI_Init(p_argc,p_argv);
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
    printf("MPI only without OpenMP\n");
  return 0;
#endif

}

void fft_init(int threads_ok)
{
  // Initialize FFTW3

#ifdef _OPENMP
  if(threads_ok)
    threads_ok=fftwf_init_threads();
  if(!threads_ok)
    msg_printf(warn,"Multi-thread FFTW not supported.\n");
#endif
    
  fftwf_mpi_init();

#ifdef _OPENMP
  if(threads_ok) {
    int nthreads=omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    msg_printf(info,"Multi-threaded FFTW: %d threads\n",nthreads);
  }
#endif

}

void snapshot_time(const float aout,const int iout, 
		   Particles const * const particles, 
		   Snapshot * const snapshot,
		   Memory mem)
{
  char filebase[256];      // TODO: make 256 to variable number...?
  timer_set_category(Snp);
  cola_set_snapshot(aout,particles,snapshot); //DAM check this
  timer_start(write);

  //Do wrap-up and relocation
  move_particles_snap(snapshot,Param.boxsize,mem.mem1,mem.size1);

  // Gadget snapshot for all particles
  // periodic wrapup not done, what about doing move_particle_min here?
  sprintf(filebase,"%s/%s%03d",Param.snap_dir,Param.snap_prefix,iout);
  write_snapshot(filebase,snapshot); //DAM: check this

  timer_stop(write);

  const double z_out=1.0/aout-1.0;
  msg_printf(normal,"snapshot %d written z=%4.2f a=%5.3f\n", 
	     iout+1,z_out,aout);
  timer_set_category(COLA);
}
