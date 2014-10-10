#ifndef _COMMON_H_
#define _COMMON_H_

#include <mpi.h>

//////
// Defined in parameters.h
typedef struct {
  //Read from param file
  int n_grid;
  double np_alloc_factor;
  int ntimestep;
  int loglevel;
  char snap_prefix[256];
  char init_prefix[256];
  char snap_dir[256];
  char init_dir[256];

  //Read from times file
  int n_aout;
  double* aout;

  //Read from IC
  unsigned long long n_part;
  double boxsize;
  double omega_m, h;
  int n_files_init;
  double a_init;
} Parameters;
Parameters Param;

int read_parameters(char *fname);


//////
// Defined in particle.h
typedef float float3[3];

typedef struct {
  float x[3];
  float dx1[3]; // ZA displacement
  float dx2[3]; // 2LPT displacement
  float v[3];   // velocity
  unsigned long long id;
} Particle;

typedef struct {
  Particle* p;
  float3* force;
  float a_x, a_v;

  int np_local, np_allocated;
  unsigned long long np_total;
  float np_average;
} Particles;

typedef struct {
  float x[3];
  float v[3];
  unsigned long long id;
} ParticleMinimum;

typedef struct {
  ParticleMinimum* p;
  int np_local;
  int np_allocated;
  unsigned long long np_total;
  float np_average;
  float a;
} Snapshot;


//////
// Defined in cola.h
void cola_kick(Particles* const particles,const float avel1);
void cola_drift(Particles* const particles,const float apos1);
void cola_set_snapshot(const double aout,Particles const * const particles,
		       Snapshot* const snapshot);

//////
// Defined in snap_io.c
typedef struct {
  int    np[6];
  double mass[6];
  double time;
  double redshift;
  int    flag_sfr;
  int    flag_feedback;
  unsigned int np_total[6];
  int    flag_cooling;
  int    num_files;
  double boxsize;
  double omega0;
  double omega_lambda;
  double hubble_param; 
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  char fill[60];
} GadgetHeader;

void print_header(GadgetHeader header);
void read_initial_condition(Particles *p);
Particle *read_snapshot_multiple(char *dirname,char *prefix,
				 unsigned long long *np_out,
				 GadgetHeader *h);
GadgetHeader read_header_multiple(char *dirname,char *prefix);
void write_snapshot(const char filebase[],Snapshot const * const snapshot);


//////
// Defined in comm.h
enum Direction {ToRight=0, ToLeft=1};

void comm_init(const int nc_pm,const float boxsize);
int comm_this_node(void);
int comm_reduce_int(int x,MPI_Op op);
int comm_share_int(int x,MPI_Op op);
float comm_xleft(const int dix);
float comm_xright(const int dix);
int comm_node(const int dix);
int comm_nnode(void);


//////
// Defined in cosmo.h
void cosmo_init(const double omega_m,const double omega_lambda);
double GrowthFactor(const double a);
double GrowthFactor2(const double a);
double Vgrowth(const double a);
double Vgrowth2(const double a);
double Qfactor(const double a);


//////
// Defined in pm.h
void pm_init(const int nc_pm,const unsigned long long np_total,
	     const float boxsize,
	     void* const mem1,const size_t size1,
	     void* const mem2,const size_t size2);
void pm_calculate_forces(Particles*);
void pm_finalize(void);


//////
// Defined in move.h
void move_particles2(Particles* const particles, const float BoxSize,
		     void* const buf, const size_t size);
void move_particles_snap(Snapshot* const particles, const float BoxSize,
			 void* const buf, const size_t size);


//////
// Defined in mem.h
typedef struct {
  void *mem1, *mem2;
  size_t size1, size2;
} Memory;

Particles* allocate_particles(const unsigned long long np_total,
			      const double np_alloc_factor);
Snapshot* allocate_snapshot(const unsigned long long np_total,const int np_alloc,
			    void* const mem,const size_t mem_size);
void allocate_shared_memory(const int nc,Memory* const mem);


//////
// Defined in msg.h
enum LogLevel {verbose,debug,normal,info,warn,error,fatal,silent};

void msg_init(void);
void msg_set_loglevel(const enum LogLevel log_level);
void msg_printf(const enum LogLevel level,const char *fmt, ...);
void msg_abort(const int errret,const char *fmt, ...);


//////
// Defined in timer.h
enum Category {Init,Reading,COLA,Snp};
enum SubCategory {all,fft,assign,force_mesh,pforce,check,
		  comm,evolve,write,kd_build,kd_link,
		  interp,global,sub};

void timer_set_category(enum Category new_cat);
void timer_start(enum SubCategory sub);
void timer_stop(enum SubCategory sub);
void timer_print();

#endif //_COMMON_H_
